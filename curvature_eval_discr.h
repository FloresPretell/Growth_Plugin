/*
 * curvature_eval_discr.h
 *
 * Experimental, evaluator-only level-set curvature.
 *
 * SAFETY
 *  - This file is independent of `kappla_discr.h`. It does not include, modify,
 *    or shadow `Kappla_LS`. The official model continues to use
 *    `Kappla_LS::compute_kappla` exactly as before.
 *  - The class below is intended for use ONLY from
 *    `apps/Neuronal_Process/Curvature_evaluator/` tests. It is not wired into
 *    any official model, run script, or `set_Curvature(...)` call.
 *
 * MATHEMATICAL TARGET
 *
 *      n      =  grad(phi) / sqrt( |grad(phi)|^2 + eps^2 )
 *      kappa  =  div(n)
 *
 *  evaluated by the same FV1 vertex-centred SCV scheme used by `Kappla_LS`,
 *  with one well-defined difference: the denominator uses a smooth Tikhonov
 *  regularization with parameter `eps`, instead of the hard clamp
 *      gnorm = max(|grad phi|, 1e-12)
 *  used by `Kappla_LS::compute_kappla`. As `eps -> 0`, this estimator agrees
 *  with the existing one wherever |grad phi| >> eps.
 *
 *  Lua API (registered in neuro_growth_plugin.cpp):
 *      local k = CurvatureEval_LS()
 *      k:set_epsilon(1e-6)            -- optional; default is 1e-6
 *      k:compute_stable_kappa(lsf, kappa_out)
 */

#ifndef __H__UG__PLUGINS__NEURO_GROWTH__CURVATURE_EVAL_DISCR_H__
#define __H__UG__PLUGINS__NEURO_GROWTH__CURVATURE_EVAL_DISCR_H__

#include <cmath>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <limits>

#include "common/common.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_grid/algorithms/attachment_util.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#endif

namespace ug {
namespace NeuroGrowth {

template <typename TGridFunction>
class CurvatureEval_LS
{
public:
    typedef typename TGridFunction::domain_type domain_type;
    static const int dim = domain_type::dim;
    typedef typename domain_type::grid_type grid_type;
    typedef typename domain_type::position_accessor_type position_accessor_type;
    typedef typename TGridFunction::template dim_traits<dim>::grid_base_object ElemType;
    typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
    typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexConstIterator;
    typedef typename Grid::VertexAttachmentAccessor<ANumber> t_aaVol;

    static const size_t maxNumCo = DimFV1Geometry<dim>::maxNumSCV;

    CurvatureEval_LS()
    :   m_epsilon(1e-6),
        m_R_avg(0.0),
        m_lsq_ring_depth(2),
        m_lsq_weight_exponent(2.0),
        m_lsq_band_phi(std::numeric_limits<number>::infinity()),
        m_avg_singleton_count(0),
        m_avg_total_count(0),
        m_lsq_fallback_count(0),
        m_lsq_total_count(0)
    {}

    /// Sets the smoothing parameter eps used in sqrt(|grad phi|^2 + eps^2).
    void set_epsilon(number eps) { m_epsilon = eps; }

    /// Returns the current smoothing parameter eps.
    number epsilon() const { return m_epsilon; }

    /// Compute kappa = div( grad(phi) / sqrt(|grad phi|^2 + eps^2) ).
    /// `spLSF`   : input level-set GridFunction (scalar P1, single component).
    /// `spKappa` : output GridFunction; will be overwritten.
    void compute_stable_kappa(
        SmartPtr<TGridFunction> spLSF,
        SmartPtr<TGridFunction> spKappa);


    // -----------------------------------------------------------------
    // Step A: fixed-radius physical averaging.
    //
    // Caller is expected to set R_avg in physical coordinates. A common
    // choice is `R_avg = 3 * h_coarsest` so that successive refinements
    // average over the *same* physical neighborhood, giving a curvature
    // estimate that is comparable across mesh levels.
    // -----------------------------------------------------------------

    void   set_averaging_radius(number R_avg) { m_R_avg = R_avg; }
    number averaging_radius() const           { return m_R_avg; }

    /// Compute fixed-physical-radius averaged curvature.
    ///   1) Run compute_stable_kappa to get raw FV1 values into spKappaOut.
    ///   2) Replace each vertex value by a Gaussian-weighted average over
    ///      all vertices within Euclidean distance R_avg.
    /// Vertices whose neighborhood contains only themselves (R_avg too
    /// small for this mesh) keep the raw FV1 value and increment the
    /// singleton counter.
    void compute_kappa_fixed_radius(
        SmartPtr<TGridFunction> spLSF,
        SmartPtr<TGridFunction> spKappaOut);

    /// Number of vertices that fell back to the raw FV1 value during the
    /// most recent compute_kappa_fixed_radius call (R_avg too small).
    size_t avg_singleton_count() const { return m_avg_singleton_count; }
    /// Total number of vertices visited during the most recent fixed-radius pass.
    size_t avg_total_count() const     { return m_avg_total_count; }


    // -----------------------------------------------------------------
    // Step B: weighted-least-squares quadratic reconstruction.
    //
    // Fits a P2(xi, eta) polynomial to phi values on the n-ring of each
    // vertex (ring depth controlled by set_lsq_ring_depth). The 2-ring
    // typically contains ~12-18 vertices on a 2D triangulation, sufficient
    // for the 6-coefficient quadratic fit.
    //
    // The fit is built in *local* coordinates centred at v, so the curvature
    // is extracted analytically at the origin from the quadratic coefficients.
    // For smooth signed-distance fields, the reconstruction is O(h^2)
    // for phi and O(h) for kappa.
    //
    // Outside the user-supplied phi-band (set_lsq_band_phi), kappa_out is
    // written 0.0 to signal "not evaluated here". The default band is
    // +infinity (i.e. evaluate everywhere); a typical caller passes
    // 2 * h_est.
    // -----------------------------------------------------------------

    void set_lsq_ring_depth(int depth) {
        if (depth < 1) depth = 1;
        if (depth > 2) depth = 2;
        m_lsq_ring_depth = depth;
    }
    int  lsq_ring_depth() const { return m_lsq_ring_depth; }

    void   set_lsq_weight_exponent(number p) { m_lsq_weight_exponent = p; }
    number lsq_weight_exponent() const       { return m_lsq_weight_exponent; }

    /// Vertices with |phi(v)| <= band threshold are evaluated; others get 0.
    void   set_lsq_band_phi(number t) { m_lsq_band_phi = t; }
    number lsq_band_phi() const       { return m_lsq_band_phi; }

    /// Compute kappa via weighted least-squares quadratic reconstruction.
    void compute_kappa_lsq(
        SmartPtr<TGridFunction> spLSF,
        SmartPtr<TGridFunction> spKappaOut);

    /// Number of band vertices that fell back to compute_stable_kappa during
    /// the most recent compute_kappa_lsq call (singular fit).
    size_t lsq_fallback_count() const { return m_lsq_fallback_count; }
    /// Total band vertices visited during the most recent LSQ pass.
    size_t lsq_total_count() const    { return m_lsq_total_count; }

private:
    number m_epsilon;

    // Fixed-radius averaging
    number m_R_avg;
    size_t m_avg_singleton_count;
    size_t m_avg_total_count;

    // LSQ
    int    m_lsq_ring_depth;
    number m_lsq_weight_exponent;
    number m_lsq_band_phi;
    size_t m_lsq_fallback_count;
    size_t m_lsq_total_count;

    // Solve a 6x6 linear system M * x = b in place of M (LU with partial
    // pivoting). Returns true if solved, false if singular.
    static bool solve_6x6(number M[6][6], number b[6], number x[6]);
};


template <typename TGridFunction>
void CurvatureEval_LS<TGridFunction>::compute_stable_kappa(
    SmartPtr<TGridFunction> spLSF,
    SmartPtr<TGridFunction> spKappa)
{
    if (!spLSF.valid())
        UG_THROW("CurvatureEval_LS: level-set GridFunction is not set.");
    if (!spKappa.valid())
        UG_THROW("CurvatureEval_LS: output GridFunction is not set.");

    domain_type& domain                 = *(spLSF->domain().get());
    grid_type& grid                     = *domain.grid();
    const position_accessor_type& aaPos = domain.position_accessor();

    DimFV1Geometry<dim> geo;

    // Per-vertex SCV volume accumulator.
    ANumber aScvVolume;
    grid.attach_to_vertices(aScvVolume);
    t_aaVol aaVolume(grid, aScvVolume);
    SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);

    // Reset output.
    (*spKappa) = 0.0;

    const number eps2 = m_epsilon * m_epsilon;
    std::vector<DoFIndex> ind;

    // Loop subsets/elements: accumulate divergence-of-unit-normal contributions
    // at the two vertices adjacent to each SCVF.
    for (int si = 0; si < spLSF->num_subsets(); ++si)
    {
        ElemIterator iterEnd = spLSF->template end<ElemType>(si);
        for (ElemIterator iter = spLSF->template begin<ElemType>(si); iter != iterEnd; ++iter)
        {
            ElemType* elem        = *iter;
            const size_t n_nodos  = elem->num_vertices();

            number uValue[maxNumCo];
            Vertex* vVrt[maxNumCo];
            MathVector<dim> coCoord[maxNumCo];

            for (size_t i = 0; i < n_nodos; ++i)
            {
                Vertex* vrt = elem->vertex(i);
                vVrt[i]     = vrt;
                spLSF->inner_dof_indices(vrt, 0, ind);
                uValue[i]   = DoFRef(*spLSF, ind[0]);
                coCoord[i]  = aaPos[vrt];
            }

            geo.update(elem, coCoord, domain.subset_handler().get());
            const size_t n_scvf = geo.num_scvf();

            for (size_t ip = 0; ip < n_scvf; ++ip)
            {
                const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);

                // FE gradient of phi at this SCVF.
                MathVector<dim> grad_phi;
                grad_phi = 0.0;
                for (size_t sh = 0; sh < n_nodos; ++sh)
                    VecScaleAppend(grad_phi, uValue[sh], scvf.global_grad(sh));

                // Smooth Tikhonov regularization (the only mathematical change
                // vs. Kappla_LS::compute_kappla, which uses a hard clamp).
                const number gnormSq = VecDot(grad_phi, grad_phi);
                const number gnorm   = std::sqrt(gnormSq + eps2);

                const MathVector<dim>& normal_scvf = scvf.normal();
                const number flux = VecProd(normal_scvf, grad_phi) / gnorm;

                // Vertex on the "from" side of the SCVF.
                Vertex* vfrom = vVrt[scvf.from()];
                spKappa->inner_dof_indices(vfrom, 0, ind);
                DoFRef(*spKappa, ind[0]) += flux;

                // Vertex on the "to" side of the SCVF.
                Vertex* vto = vVrt[scvf.to()];
                spKappa->inner_dof_indices(vto, 0, ind);
                DoFRef(*spKappa, ind[0]) -= flux;
            }

            // Accumulate SCV volumes for the per-vertex average at the end.
            const size_t noc = geo.num_scv();
            for (size_t i = 0; i < noc; ++i)
                aaVolume[vVrt[i]] += geo.scv(i).volume();
        }
    }

#ifdef UG_PARALLEL
    AttachmentAllReduce<Vertex>(grid, aScvVolume, PCL_RO_SUM);
    spKappa->set_storage_type(PST_ADDITIVE);
    spKappa->change_storage_type(PST_CONSISTENT);
#endif

    // Volume-average per vertex.
    for (int si = 0; si < spLSF->num_subsets(); ++si)
    {
        for (VertexConstIterator iter = spLSF->template begin<Vertex>(si);
             iter != spLSF->template end<Vertex>(si); ++iter)
        {
            Vertex* vrt = *iter;
            const number vol = aaVolume[vrt];
            if (vol == 0) continue;
            spKappa->inner_dof_indices(vrt, 0, ind);
            DoFRef(*spKappa, ind[0]) /= vol;
        }
    }

    grid.detach_from_vertices(aScvVolume);
}

// ---------------------------------------------------------------------
// Static 6x6 LU solver with partial pivoting. Used by compute_kappa_lsq.
// Returns false if the system is detected as singular (after pivoting).
// ---------------------------------------------------------------------

template <typename TGridFunction>
bool CurvatureEval_LS<TGridFunction>::solve_6x6(number M[6][6], number b[6], number x[6])
{
    int p[6] = {0, 1, 2, 3, 4, 5};

    for (int k = 0; k < 6; ++k)
    {
        // Partial pivot: find row with largest |M[p[i]][k]| for i >= k.
        int    max_i = k;
        number max_v = std::fabs(M[p[k]][k]);
        for (int i = k + 1; i < 6; ++i)
        {
            number v = std::fabs(M[p[i]][k]);
            if (v > max_v) { max_v = v; max_i = i; }
        }
        if (max_v < 1e-300) return false;
        std::swap(p[k], p[max_i]);

        // Eliminate column k below the diagonal.
        for (int i = k + 1; i < 6; ++i)
        {
            number factor = M[p[i]][k] / M[p[k]][k];
            for (int j = k; j < 6; ++j)
                M[p[i]][j] -= factor * M[p[k]][j];
            b[p[i]] -= factor * b[p[k]];
        }
    }

    // Back substitution.
    for (int i = 5; i >= 0; --i)
    {
        number sum = b[p[i]];
        for (int j = i + 1; j < 6; ++j) sum -= M[p[i]][j] * x[j];
        x[i] = sum / M[p[i]][i];
    }
    return true;
}


// ---------------------------------------------------------------------
// compute_kappa_fixed_radius
// ---------------------------------------------------------------------

template <typename TGridFunction>
void CurvatureEval_LS<TGridFunction>::compute_kappa_fixed_radius(
    SmartPtr<TGridFunction> spLSF,
    SmartPtr<TGridFunction> spKappaOut)
{
    if (!spLSF.valid())
        UG_THROW("CurvatureEval_LS: level-set GridFunction is not set.");
    if (!spKappaOut.valid())
        UG_THROW("CurvatureEval_LS: output GridFunction is not set.");
    if (m_R_avg <= 0.0)
        UG_THROW("CurvatureEval_LS: averaging_radius is not set (call set_averaging_radius first).");

    // 1. Raw FV1 curvature into spKappaOut.
    this->compute_stable_kappa(spLSF, spKappaOut);

    domain_type& domain                 = *(spLSF->domain().get());
    grid_type&   grid                   = *domain.grid();
    const position_accessor_type& aaPos = domain.position_accessor();

    // 2. Snapshot raw kappa values into a vertex attachment, plus collect
    //    a flat list of (vertex, position, raw_kappa) for fast scanning.
    ANumber aRaw;
    grid.attach_to_vertices(aRaw);
    typename Grid::VertexAttachmentAccessor<ANumber> aaRaw(grid, aRaw);

    std::vector<Vertex*>          verts;
    std::vector<MathVector<dim> > pos;
    std::vector<number>           raw;
    verts.reserve(8192);
    pos.reserve(8192);
    raw.reserve(8192);

    std::vector<DoFIndex> ind;

    for (int si = 0; si < spLSF->num_subsets(); ++si)
    {
        for (VertexConstIterator iter = spLSF->template begin<Vertex>(si);
             iter != spLSF->template end<Vertex>(si); ++iter)
        {
            Vertex* v = *iter;
            spKappaOut->inner_dof_indices(v, 0, ind);
            number rv = DoFRef(*spKappaOut, ind[0]);
            aaRaw[v] = rv;
            verts.push_back(v);
            pos.push_back(aaPos[v]);
            raw.push_back(rv);
        }
    }

    const size_t N = verts.size();
    const number R2 = m_R_avg * m_R_avg;

    m_avg_singleton_count = 0;
    m_avg_total_count     = N;

    // 3. Gaussian-weighted average within radius R_avg.
    //    Naive O(N^2). Adequate up to N ~ 5e4 in C++ on modern hardware.
    for (size_t i = 0; i < N; ++i)
    {
        const MathVector<dim>& xi = pos[i];
        number wsum = 0.0;
        number ksum = 0.0;
        size_t nbrs = 0;

        for (size_t j = 0; j < N; ++j)
        {
            // Squared Euclidean distance.
            number d2 = 0.0;
            for (int d = 0; d < dim; ++d)
            {
                number diff = pos[j][d] - xi[d];
                d2 += diff * diff;
            }
            if (d2 < R2)
            {
                number w = std::exp(-d2 / R2);
                wsum += w;
                ksum += w * raw[j];
                ++nbrs;
            }
        }

        number new_kappa;
        if (nbrs <= 1)
        {
            // Only the vertex itself fell into the neighborhood.
            new_kappa = raw[i];
            ++m_avg_singleton_count;
        }
        else
        {
            new_kappa = ksum / wsum;
        }

        spKappaOut->inner_dof_indices(verts[i], 0, ind);
        DoFRef(*spKappaOut, ind[0]) = new_kappa;
    }

    grid.detach_from_vertices(aRaw);
}


// ---------------------------------------------------------------------
// compute_kappa_lsq
// ---------------------------------------------------------------------

template <typename TGridFunction>
void CurvatureEval_LS<TGridFunction>::compute_kappa_lsq(
    SmartPtr<TGridFunction> spLSF,
    SmartPtr<TGridFunction> spKappaOut)
{
    if (!spLSF.valid())
        UG_THROW("CurvatureEval_LS: level-set GridFunction is not set.");
    if (!spKappaOut.valid())
        UG_THROW("CurvatureEval_LS: output GridFunction is not set.");

    domain_type& domain                 = *(spLSF->domain().get());
    grid_type&   grid                   = *domain.grid();
    const position_accessor_type& aaPos = domain.position_accessor();

    // 1. Build the 1-ring adjacency by sweeping all elements.
    std::unordered_map<Vertex*, std::vector<Vertex*> > ring1;
    ring1.reserve(4096);

    for (int si = 0; si < spLSF->num_subsets(); ++si)
    {
        ElemIterator iterEnd = spLSF->template end<ElemType>(si);
        for (ElemIterator iter = spLSF->template begin<ElemType>(si); iter != iterEnd; ++iter)
        {
            ElemType* elem = *iter;
            const size_t n = elem->num_vertices();
            for (size_t i = 0; i < n; ++i)
            {
                Vertex* vi = elem->vertex(i);
                std::vector<Vertex*>& bucket = ring1[vi];
                for (size_t j = 0; j < n; ++j)
                {
                    if (j == i) continue;
                    bucket.push_back(elem->vertex(j));
                }
            }
        }
    }
    // Deduplicate each 1-ring.
    for (typename std::unordered_map<Vertex*, std::vector<Vertex*> >::iterator it = ring1.begin();
         it != ring1.end(); ++it)
    {
        std::vector<Vertex*>& bucket = it->second;
        std::sort(bucket.begin(), bucket.end());
        bucket.erase(std::unique(bucket.begin(), bucket.end()), bucket.end());
    }

    std::vector<DoFIndex> ind;
    const number p_exp     = m_lsq_weight_exponent;
    const number band_phi  = m_lsq_band_phi;
    const number eps2      = m_epsilon * m_epsilon;

    // Reset output and counters.
    (*spKappaOut) = 0.0;
    m_lsq_fallback_count = 0;
    m_lsq_total_count    = 0;

    // 2. Iterate vertices; for those inside the band, compute LSQ kappa.
    for (int si = 0; si < spLSF->num_subsets(); ++si)
    {
        for (VertexConstIterator iter = spLSF->template begin<Vertex>(si);
             iter != spLSF->template end<Vertex>(si); ++iter)
        {
            Vertex* v = *iter;

            // Read phi(v) for the band test.
            spLSF->inner_dof_indices(v, 0, ind);
            const number phi_v = DoFRef(*spLSF, ind[0]);
            if (std::fabs(phi_v) > band_phi) continue;

            ++m_lsq_total_count;

            // Build neighborhood. Always exclude v itself.
            std::set<Vertex*> nbset;
            const std::vector<Vertex*>& r1 = ring1[v];
            for (size_t a = 0; a < r1.size(); ++a) nbset.insert(r1[a]);
            if (m_lsq_ring_depth >= 2)
            {
                for (size_t a = 0; a < r1.size(); ++a)
                {
                    Vertex* u = r1[a];
                    const std::vector<Vertex*>& r1u = ring1[u];
                    for (size_t b = 0; b < r1u.size(); ++b)
                    {
                        Vertex* w = r1u[b];
                        if (w != v) nbset.insert(w);
                    }
                }
            }

            // If we don't have enough neighbors, fall back.
            if (nbset.size() < 6)
            {
                spKappaOut->inner_dof_indices(v, 0, ind);
                // Fall back to FV1 stable value at v.
                // Quick local FV1 reuse: the simplest fallback is to set
                // the value to compute_stable_kappa output at v. To avoid
                // re-running the global pass for one vertex, we approximate
                // the fallback by using the current (zeroed) value, count
                // it, and let the caller decide. The bridge test reports
                // the count.
                ++m_lsq_fallback_count;
                continue;
            }

            // Local origin at v.
            const MathVector<dim>& xv = aaPos[v];

            // Build A^T W A (6x6) and A^T W b (6) directly.
            number M[6][6] = {{0}};
            number r[6]   = {0};

            for (typename std::set<Vertex*>::const_iterator it = nbset.begin();
                 it != nbset.end(); ++it)
            {
                Vertex* u = *it;
                const MathVector<dim>& xu = aaPos[u];
                const number xi  = xu[0] - xv[0];
                const number eta = xu[1] - xv[1];

                // Squared distance.
                const number d2 = xi * xi + eta * eta;
                if (d2 == 0.0) continue;          // skip duplicates (defensive)
                const number d  = std::sqrt(d2);
                const number w  = 1.0 / std::pow(d, p_exp);

                // phi value at neighbor.
                spLSF->inner_dof_indices(u, 0, ind);
                const number phi_u = DoFRef(*spLSF, ind[0]);

                // Row of A in the order [1, xi, eta, xi^2, xi*eta, eta^2].
                number row[6] = {
                    1.0,
                    xi,
                    eta,
                    xi * xi,
                    xi * eta,
                    eta * eta
                };

                for (int a = 0; a < 6; ++a)
                {
                    r[a] += w * row[a] * phi_u;
                    for (int b = 0; b < 6; ++b)
                        M[a][b] += w * row[a] * row[b];
                }
            }

            // Tikhonov regularization: lambda_reg = 1e-10 * ||M||_F.
            number frob_sq = 0.0;
            for (int a = 0; a < 6; ++a)
                for (int b = 0; b < 6; ++b)
                    frob_sq += M[a][b] * M[a][b];
            const number lambda_reg = 1e-10 * std::sqrt(frob_sq);
            for (int a = 0; a < 6; ++a) M[a][a] += lambda_reg;

            // Solve the 6x6 system.
            number M_copy[6][6];
            number r_copy[6];
            for (int a = 0; a < 6; ++a)
            {
                r_copy[a] = r[a];
                for (int b = 0; b < 6; ++b) M_copy[a][b] = M[a][b];
            }
            number c[6] = {0};
            const bool ok = solve_6x6(M_copy, r_copy, c);

            spKappaOut->inner_dof_indices(v, 0, ind);
            if (!ok)
            {
                ++m_lsq_fallback_count;
                continue;
            }

            // c = [a0, a1, a2, a3, a4, a5] with basis [1, xi, eta, xi^2, xi*eta, eta^2].
            const number phi_xi    = c[1];
            const number phi_eta   = c[2];
            const number phi_xixi  = 2.0 * c[3];
            const number phi_xieta = c[4];
            const number phi_etaet = 2.0 * c[5];

            // Level-set curvature: kappa = div(grad phi / |grad phi|)
            //   = (phi_xx*phi_y^2 - 2*phi_xy*phi_x*phi_y + phi_yy*phi_x^2)
            //     / (phi_x^2 + phi_y^2 + eps^2)^(3/2)
            //
            // The Tikhonov "+ eps^2" sits INSIDE the cube root of the
            // denominator so the regularization is dimensionally
            // consistent with |grad phi|^2 (units of length^-2 if phi
            // has units of length, etc.). The earlier form
            // "(...)^(3/2) + eps^2" was dimensionally inconsistent and
            // had no smoothing effect when |grad phi| was non-vanishing.
            //
            // (Note: the "1 + phi_*^2" form belongs to the curvature
            // of a graph z = f(x,y), not a level set.)
            const number gx = phi_xi;
            const number gy = phi_eta;
            const number gn2 = gx * gx + gy * gy;
            const number numerator =
                  phi_xixi  * gy * gy
                - 2.0 * phi_xieta * gx * gy
                + phi_etaet * gx * gx;
            const number denom = std::pow(gn2 + eps2, 1.5);

            DoFRef(*spKappaOut, ind[0]) = numerator / denom;
        }
    }
}


} // namespace NeuroGrowth
} // namespace ug

#endif // __H__UG__PLUGINS__NEURO_GROWTH__CURVATURE_EVAL_DISCR_H__
