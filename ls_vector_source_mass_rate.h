/*
 * ls_vector_source_mass_rate.h
 *
 * Diagnostic helper: computes the net mass rate contributed by a vector
 * source S to the intracellular subdomain (lsf <= 0) using the same FV1
 * discrete convention as ConvectionDiffusionFV1::add_def_A_elem:
 *
 *   flux_ip  = VecDot(S[ip], scvf.normal())
 *   d[from] -= flux_ip
 *   d[to]   += flux_ip
 *
 * The net mass rate into Omega_i (lsf <= 0) is the sum of contributions
 * to interior nodes only. Internal SCVFs (both nodes in same domain) cancel.
 * Only interface SCVFs (one node on each side) contribute:
 *
 *   net_rate += -flux  if from-node in Omega_i, to-node in Omega_e
 *   net_rate +=  flux  if from-node in Omega_e, to-node in Omega_i
 *
 * Intended use: compute synthesis_U_rate_direct from LSInfluxTopBranch.
 */

#ifndef __H__UG__PLUGINS__LS_VECTOR_SOURCE_MASS_RATE_H__
#define __H__UG__PLUGINS__LS_VECTOR_SOURCE_MASS_RATE_H__

#include <vector>

#include "common/common.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/domain_util.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif

namespace ug {
namespace NeuroGrowth {

/**
 * FVLSVectorSourceIntegral
 *
 * Computes the net FV1 mass rate from a vector CplUserData into lsf <= 0.
 *
 * Lua API (example):
 *   local synth_integrator = FVLSVectorSourceIntegral(gf.lsf)
 *   synth_integrator:set_vector_source(Syntesis_U_linker)
 *   synth_integrator:compute(time)
 *   local synthesis_rate = synth_integrator:result()
 */
template <typename TGridFunction>
class FVLSVectorSourceIntegral
{
public:
    typedef TGridFunction gf_type;
    static const int dim = gf_type::dim;
    typedef typename gf_type::domain_type domain_type;
    typedef typename domain_type::position_accessor_type position_accessor_type;

    // LSF grid function type: always CPUAlgebra (scalar)
    typedef CPUAlgebra lsf_algebra_type;
    typedef GridFunction<domain_type, lsf_algebra_type> ls_gf_type;

    // Vector source type
    typedef CplUserData<MathVector<dim>, dim> vec_src_type;

public:
    FVLSVectorSourceIntegral(SmartPtr<ls_gf_type> sp_lsf)
    :   m_spLSF(sp_lsf), m_time(0), m_result(0)
    {}

    virtual ~FVLSVectorSourceIntegral() {}

    void set_vector_source(SmartPtr<vec_src_type> sp_src)
    {
        m_spVecSrc = sp_src;
    }

    /// Compute net mass rate. Call once per time step before reading result().
    void compute(number time)
    {
        if (!m_spVecSrc.valid())
            UG_THROW("FVLSVectorSourceIntegral: No vector source set. "
                     "Call set_vector_source() first.");
        if (!m_spLSF.valid())
            UG_THROW("FVLSVectorSourceIntegral: No LSF grid function.");

        typedef typename domain_traits<dim>::DimElemList ElemList;
        m_time   = time;
        m_result = 0.0;
        boost::mpl::for_each<ElemList>(AddContributions(this));

#ifdef UG_PARALLEL
        pcl::ProcessCommunicator procComm;
        m_result = procComm.allreduce(m_result, PCL_RO_SUM);
#endif
    }

    /// Returns the net mass rate into Omega_i (lsf <= 0) after compute().
    number result() const { return m_result; }

private:
    template <typename TElem>
    void add_contributions_of_all()
    {
        typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_t;
        typedef FV1Geometry<TElem, dim> TFVGeom;
        static const size_t num_corners = ref_elem_t::numCorners;

        const ls_gf_type& lsf   = *m_spLSF;
        const position_accessor_type& aaPos = lsf.domain()->position_accessor();

        std::vector<DoFIndex> ind(1);
        MathVector<dim>  corners[num_corners];
        number           lsf_values[num_corners];

        for (int si = 0; si < lsf.num_subsets(); si++)
        {
            typedef typename ls_gf_type::template traits<TElem>::const_iterator ElemIter;
            ElemIter iterEnd = lsf.template end<TElem>(si);

            for (ElemIter iter = lsf.template begin<TElem>(si); iter != iterEnd; ++iter)
            {
                TElem* elem = *iter;

                // Collect corner coordinates and LSF nodal values
                bool has_neg = false, has_pos = false;
                for (size_t i = 0; i < num_corners; i++)
                {
                    Vertex* vrt = elem->vertex(i);
                    corners[i]  = aaPos[vrt];

                    if (lsf.inner_dof_indices(vrt, 0, ind) != 1)
                        UG_THROW("FVLSVectorSourceIntegral: LSF must be a scalar P1 function.");

                    lsf_values[i] = DoFRef(lsf, ind[0]);
                    if (lsf_values[i] <= 0.0) has_neg = true;
                    else                      has_pos = true;
                }

                // Skip elements fully inside or outside (not cut)
                if (!has_neg || !has_pos) continue;

                // Update FV1 geometry for this cut element
                static TFVGeom& geo = GeomProvider<TFVGeom>::get();
                try { geo.update(elem, corners); }
                UG_CATCH_THROW("FVLSVectorSourceIntegral: Cannot update FV geometry.");

                const size_t num_scvf = geo.num_scvf();

                // Evaluate vector source at SCVF midpoints
                // Uses the same operator() calling convention as Export.h (line 1757)
                std::vector<MathVector<dim>> vVal(num_scvf);
                (*m_spVecSrc)(vVal.data(),
                              geo.scvf_global_ips(),
                              m_time, si,
                              elem, geo.corners(),
                              geo.scvf_local_ips(),
                              num_scvf, NULL);

                // Accumulate net flux contribution to Omega_i (lsf <= 0)
                for (size_t ip = 0; ip < num_scvf; ip++)
                {
                    const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
                    const size_t from = scvf.from();
                    const size_t to   = scvf.to();
                    const bool from_neg = (lsf_values[from] <= 0.0);
                    const bool to_neg   = (lsf_values[to]   <= 0.0);

                    // Internal SCVF: both nodes in same domain → cancels for mass total
                    if (from_neg == to_neg) continue;

                    const number flux = VecDot(vVal[ip], scvf.normal());

                    // d[from] -= flux  →  Omega_i[from] gets -flux
                    // d[to]   += flux  →  Omega_i[to]   gets +flux
                    if (from_neg) m_result -= flux;
                    else          m_result += flux;
                }
            }
        }
    }

    // Functor for boost::mpl::for_each over element types
    struct AddContributions
    {
        FVLSVectorSourceIntegral* m_pThis;
        explicit AddContributions(FVLSVectorSourceIntegral* p) : m_pThis(p) {}
        template <typename TElem>
        void operator()(TElem) { m_pThis->template add_contributions_of_all<TElem>(); }
    };

private:
    SmartPtr<ls_gf_type>   m_spLSF;
    SmartPtr<vec_src_type> m_spVecSrc;
    number                 m_time;
    number                 m_result;
};

} // namespace NeuroGrowth
} // namespace ug

#endif // __H__UG__PLUGINS__LS_VECTOR_SOURCE_MASS_RATE_H__
