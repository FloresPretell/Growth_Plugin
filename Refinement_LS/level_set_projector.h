/*
 * LevelSetProjector — a refinement projector guided by a level-set function φ.
 *
 * Generalises SphereProjector from an implicit SPHERE (φ = |x−c| − R) to an
 * arbitrary implicit surface {φ = 0}. When refinement creates a new vertex, the
 * projector places it on the isolevel through the parent-average level (the exact
 * analogue of SphereProjector projecting onto the sphere through the parent-average
 * radius). For a vertex born on a cut edge (parents straddle φ=0) the isolevel is
 * φ=0, so the reconstructed interface stays on the TRUE surface, not the coarse
 * polygon.
 *
 * v1: ANALYTIC φ (sphere / cylinder), evaluated through the general Newton-on-
 * isolevel machinery. Reproducing SphereProjector for the sphere case is the
 * correctness test. A discrete (GridFunction-attachment) φ is the planned v2.
 *
 * Works in 2D and 3D (coords are 3D internally), with any refiner (hanging or
 * conforming — new_vertex is called by all), and any mesh (position-only).
 */
#ifndef __H__UG_level_set_projector
#define __H__UG_level_set_projector

#include <functional>
#include <cmath>
#include "lib_grid/refinement/projectors/refinement_projector.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/multi_grid.h"

namespace ug {
namespace LevelSet {

class LevelSetProjector : public RefinementProjector
{
	public:
		LevelSetProjector() : m_band(-1), m_newtonIters(3), m_fdH(1e-6) {}

		LevelSetProjector(SPIGeometry3d geometry) :
			RefinementProjector(geometry), m_band(-1), m_newtonIters(3), m_fdH(1e-6) {}

	///	analytic φ = |x − c| − R  (implicit sphere/circle — reproduces SphereProjector)
		void set_sphere(const vector3& c, number R)
		{
			m_phi = [c, R](const vector3& x) -> number {
				vector3 d; VecSubtract(d, x, c); return VecLength(d) - R;
			};
		}

	///	analytic φ = dist(x, axis) − R  (implicit cylinder; axis point p, unit dir a)
		void set_cylinder(const vector3& p, const vector3& a, number R)
		{
			vector3 au = a; VecNormalize(au, au);
			m_phi = [p, au, R](const vector3& x) -> number {
				vector3 d; VecSubtract(d, x, p);
				number t = VecDot(d, au);
				vector3 ax; VecScale(ax, au, t);
				vector3 r; VecSubtract(r, d, ax);
				return VecLength(r) - R;
			};
		}

	///	only project vertices with |φ_parent| ≤ band; leave the far field linear (<0 = always)
		void set_band(number b)          { m_band = b; }
		void set_newton_iterations(int n){ m_newtonIters = n; }

	// ---- SHELL mode: organize the whole INTERIOR onto nested φ-isosurfaces --------
	///	When shell mode is on, a new interior vertex (φ_parent ≤ band) is projected
	///	onto the NEAREST discrete shell φ = k·dPhi (k integer) along −∇φ, instead of
	///	being left on the box lattice. This turns a Cartesian octree of the interior
	///	into concentric radial shells following the tube (the level-set analogue of
	///	the NeuriteProjector's fractional radial coordinate). Exterior nodes (φ>band)
	///	stay on the clean grid. dPhi = shell spacing in φ-space; 0 disables.
		void set_shell_spacing(number dPhi) { m_dPhi = dPhi; m_shellMode = (dPhi > 0); }
	///	nShells fractional-radius shells across a tube of radius R (dPhi = R/nShells)
		void set_num_shells(int nShells, number R)
		{	set_shell_spacing(nShells > 0 ? R / (number)nShells : 0.0); }
	///	inversion guard: never move a vertex farther than frac·(local spacing) from
	///	its linear position (caps radial displacement so thin cells cannot invert)
		void set_max_radial_disp(number frac) { m_maxRadFrac = frac; }

	// ---- v2: DISCRETE φ from a GridFunction, carried as grid attachments -------
	///	attach the per-vertex φ / ∇φ storage to the grid (called by the templated
	///	set-LSF helper before it fills the values). Switches to discrete mode.
		void attach_grid(Grid& g)
		{
			m_pGrid = &g;
			if (!g.has_vertex_attachment(m_aLSF))  g.attach_to_vertices(m_aLSF);
			if (!g.has_vertex_attachment(m_aGrad)) g.attach_to_vertices(m_aGrad);
			m_aaLSF.access(g, m_aLSF);
			m_aaGrad.access(g, m_aGrad);
			m_discrete = true;
		}
		void set_vertex_data(Vertex* v, number phi, const vector3& grad)
		{	m_aaLSF[v] = phi;  m_aaGrad[v] = grad; }
		bool discrete() const { return m_discrete; }

		virtual number new_vertex(Vertex* vrt, Vertex* parent)
		{	set_pos(vrt, pos(parent));
			if (m_discrete) set_vertex_data(vrt, m_aaLSF[parent], m_aaGrad[parent]);
			return 1; }

		virtual number new_vertex(Vertex* vrt, Edge* parent)
		{	return m_discrete ? project_edge_discrete(vrt, parent) : project(vrt, parent); }
		virtual number new_vertex(Vertex* vrt, Face* parent)
		{	return m_discrete ? project_center_discrete(vrt, parent) : project(vrt, parent); }
		virtual number new_vertex(Vertex* vrt, Volume* parent)
		{	return m_discrete ? project_center_discrete(vrt, parent) : project(vrt, parent); }

	private:
		number phi(const vector3& x) const { return m_phi ? m_phi(x) : 0.0; }

	///	snap a parent level φ0 to the nearest interior shell φ = k·dPhi; never past the
	///	surface into the exterior (returns 0 if the nearest shell would be positive)
		number snap_shell(number phi0) const
		{
			if (!m_shellMode || m_dPhi <= 0) return phi0;
			const number s = std::round(phi0 / m_dPhi) * m_dPhi;
			return (s > 0.0) ? 0.0 : s;
		}

		vector3 grad_phi(const vector3& x) const
		{
			vector3 g(0, 0, 0);
			for (int d = 0; d < 3; ++d) {
				vector3 xp = x, xm = x;
				xp[d] += m_fdH; xm[d] -= m_fdH;
				g[d] = (phi(xp) - phi(xm)) / (2.0 * m_fdH);
			}
			return g;
		}

		template <class TElem>
		number project(Vertex* vrt, TElem* parent)
		{
			typename TElem::ConstVertexArray vrts = parent->vertices();
			const size_t n = parent->num_vertices();
			if (n == 0) { set_pos(vrt, vector3(0, 0, 0)); return 1; }

		//	parent-average position x0 and level φ0 (SphereProjector analogue)
			vector3 x0(0, 0, 0);
			number phi0 = 0;
			for (size_t i = 0; i < n; ++i) {
				vector3 p = pos(vrts[i]);
				VecAdd(x0, x0, p);
				phi0 += phi(p);
			}
			VecScale(x0, x0, 1.0 / (number)n);
			phi0 /= (number)n;

		//	gate: shell mode projects the WHOLE interior (φ0≤band), snapping to the
		//	nearest shell; legacy mode projects the symmetric band onto {φ=φ0}. The
		//	exterior (φ0>band) always keeps the linear midpoint (clean outer grid).
			const bool skip = m_shellMode ? (m_band > 0 && phi0 > m_band)
			                              : (m_band > 0 && std::fabs(phi0) > m_band);
			if (skip) { set_pos(vrt, x0); return 1; }
			const number target = snap_shell(phi0);   // == phi0 when shell mode off

		//	Newton-project x0 onto the isolevel {φ = target} along ∇φ
			vector3 x = x0;
			for (int it = 0; it < m_newtonIters; ++it) {
				number f = phi(x) - target;
				vector3 g = grad_phi(x);
				number gg = VecDot(g, g);
				if (gg < 1e-20) break;
				vector3 step; VecScale(step, g, f / gg);
				VecSubtract(x, x, step);
			}
		//	inversion guard: cap the radial displacement to a fraction of the parent
		//	extent, else a thin cell can flip (negative Jacobian). Revert to linear.
			if (m_maxRadFrac > 0) {
				number localLen = 0;
				for (size_t i = 0; i < n; ++i) {
					const number d = VecDistance(pos(vrts[i]), x0);
					if (d > localLen) localLen = d;
				}
				if (VecDistance(x, x0) > m_maxRadFrac * localLen) x = x0;
			}
			set_pos(vrt, x);
			return 1;
		}

	// ---- discrete edge projection: cubic-Hermite φ=0 crossing along the edge ----
	//	A cut edge carries a CURVED interface (from the two endpoint φ AND ∇φ), so
	//	the new vertex lands on the true surface, not the linear midpoint. This is
	//	the reconstruction that beats a carried linear φ (the Part-C gap).
		number project_edge_discrete(Vertex* vrt, Edge* parent)
		{
			Vertex* va = parent->vertex(0);
			Vertex* vb = parent->vertex(1);
			const vector3 pa = pos(va), pb = pos(vb);
			const number  fa = m_aaLSF[va], fb = m_aaLSF[vb];

			vector3 e; VecSubtract(e, pb, pa);               // edge vector a->b
			const number ma = VecDot(m_aaGrad[va], e);       // dφ/ds at a (Hermite slope)
			const number mb = VecDot(m_aaGrad[vb], e);       // dφ/ds at b

			number s = 0.5;
			bool cut = (fa < 0.0) != (fb < 0.0);
			if (cut) {
			//	linear guess for the φ=0 crossing — ALWAYS strictly interior for a cut
			//	edge (fa,fb opposite signs), so it never coincides with an endpoint
				const number s0 = (fa != fb) ? fa / (fa - fb) : 0.5;
				s = s0;
			//	Newton on the Hermite cubic to place the vertex on the CURVED interface
				for (int it = 0; it < (m_newtonIters > 0 ? m_newtonIters : 3); ++it) {
					const number f  = hermite(s, fa, ma, fb, mb);
					const number df = hermite_d(s, fa, ma, fb, mb);
					if (std::fabs(df) < 1e-20) break;
					s -= f / df;
				}
			//	In high-curvature regions the cubic Newton can overshoot the edge;
			//	falling back to the (always-interior) linear crossing avoids placing
			//	the new vertex ON an endpoint (a zero-area sliver). Only trust the
			//	Hermite result if it stayed near the linear crossing and interior.
				if (!(s > 0.0 && s < 1.0) || std::fabs(s - s0) > 0.35)
					s = s0;
				const number smin = 0.05, smax = 0.95;   // never coincide with an endpoint
				if (s < smin) s = smin; else if (s > smax) s = smax;
			}
			else if (m_shellMode) {
			//	NON-cut interior edge: instead of the midpoint, place the new vertex on
			//	the nearest φ-shell so the interior organizes into concentric shells.
			//	Only for interior edges (φ0≤band); exterior edges keep the midpoint so
			//	the outer box grid stays regular.
				const number phi0 = 0.5 * (fa + fb);
				if (!(m_band > 0 && phi0 > m_band)) {
					const number target = snap_shell(phi0);
					const number s0 = 0.5;
					s = s0;
					for (int it = 0; it < (m_newtonIters > 0 ? m_newtonIters : 3); ++it) {
						const number f  = hermite(s, fa, ma, fb, mb) - target;
						const number df = hermite_d(s, fa, ma, fb, mb);
						if (std::fabs(df) < 1e-20) break;
						s -= f / df;
					}
				//	inversion guard: displacement from the midpoint is |s−0.5| of the
				//	edge length; cap it and revert to the midpoint if it overshoots
					if (!(s > 0.0 && s < 1.0) || std::fabs(s - s0) > m_maxRadFrac)
						s = s0;
					const number smin = 0.05, smax = 0.95;
					if (s < smin) s = smin; else if (s > smax) s = smax;
				}
			}
			vector3 x; VecScaleAdd(x, 1.0 - s, pa, s, pb);
			set_pos(vrt, x);

		//	carry φ / ∇φ onto the new vertex (0 exactly on the interface) so multi-
		//	level refinement keeps working
			const number fnew = cut ? 0.0 : hermite(s, fa, ma, fb, mb);
			vector3 gnew; VecScaleAdd(gnew, 1.0 - s, m_aaGrad[va], s, m_aaGrad[vb]);
			set_vertex_data(vrt, fnew, gnew);
			return 1;
		}

	//	face/volume centre: linear interpolate φ/∇φ from corners; one Newton step
	//	onto φ=0 if the centre is within the band (rarely on the interface for
	//	simplicial refinement, where new vertices are on edges).
		template <class TElem>
		number project_center_discrete(Vertex* vrt, TElem* parent)
		{
			typename TElem::ConstVertexArray vrts = parent->vertices();
			const size_t n = parent->num_vertices();
			if (n == 0) { set_pos(vrt, vector3(0, 0, 0)); return 1; }
			vector3 x0(0, 0, 0), g0(0, 0, 0); number f0 = 0;
			for (size_t i = 0; i < n; ++i) {
				VecAdd(x0, x0, pos(vrts[i]));
				f0 += m_aaLSF[vrts[i]];
				VecAdd(g0, g0, m_aaGrad[vrts[i]]);
			}
			VecScale(x0, x0, 1.0 / (number)n);
			f0 /= (number)n; VecScale(g0, g0, 1.0 / (number)n);

			vector3 x = x0;
		//	shell mode: march the centre toward the nearest shell (φ=target) over the
		//	whole interior; legacy: march onto φ=0 within the band. Exterior untouched.
			const bool doProj = m_shellMode ? (!(m_band > 0 && f0 > m_band))
			                                : (m_band <= 0 || std::fabs(f0) <= m_band);
			if (doProj) {
				const number target = m_shellMode ? snap_shell(f0) : 0.0;
				const number gg = VecDot(g0, g0);
				if (gg > 1e-20) {
					vector3 step; VecScale(step, g0, (f0 - target) / gg);
					VecSubtract(x, x, step);
				}
			}
		//	inversion guard: cap displacement to a fraction of the parent extent
			if (m_maxRadFrac > 0) {
				number localLen = 0;
				for (size_t i = 0; i < n; ++i) {
					const number d = VecDistance(pos(vrts[i]), x0);
					if (d > localLen) localLen = d;
				}
				if (VecDistance(x, x0) > m_maxRadFrac * localLen) x = x0;
			}
			set_pos(vrt, x);
			set_vertex_data(vrt, f0, g0);
			return 1;
		}

		static number hermite(number s, number f0, number m0, number f1, number m1)
		{
			const number s2 = s * s, s3 = s2 * s;
			return (2*s3 - 3*s2 + 1) * f0 + (s3 - 2*s2 + s) * m0
			     + (-2*s3 + 3*s2)    * f1 + (s3 - s2)       * m1;
		}
		static number hermite_d(number s, number f0, number m0, number f1, number m1)
		{
			const number s2 = s * s;
			return (6*s2 - 6*s) * f0 + (3*s2 - 4*s + 1) * m0
			     + (-6*s2 + 6*s) * f1 + (3*s2 - 2*s)     * m1;
		}

		std::function<number(const vector3&)> m_phi;
		number m_band;
		int    m_newtonIters;
		number m_fdH;

	//	shell mode (project interior onto nested φ = k·dPhi shells)
		number m_dPhi       = 0.0;
		bool   m_shellMode  = false;
		number m_maxRadFrac = 0.5;

	//	discrete (GridFunction) mode
		Grid*                                          m_pGrid = 0;
		Attachment<number>                             m_aLSF;
		Attachment<vector3>                            m_aGrad;
		Grid::VertexAttachmentAccessor<Attachment<number> >  m_aaLSF;
		Grid::VertexAttachmentAccessor<Attachment<vector3> > m_aaGrad;
		bool                                           m_discrete = false;
};

}  // namespace LevelSet
}  // namespace ug

#endif
