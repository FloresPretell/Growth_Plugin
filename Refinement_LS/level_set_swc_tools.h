/*
 * InterpolateSWCPhi — fill a level-set GridFunction with φ = min over SWC segments of
 * (dist(x, segment) − radius − radiusOffset), in C++ (no per-vertex Lua loop). This is
 * the fast replacement for the Lua `Interpolate("PhiFld", u, ...)` in the growth loop:
 * moving the SWC min-distance evaluation into C++ makes the adaptive growth benefit
 * from the Release build (the Lua interpreter was the bottleneck). `radiusOffset` grows
 * the tube (inflation): larger offset ⇒ larger radius ⇒ interface advances outward.
 */
#ifndef __H__UG_level_set_swc_tools
#define __H__UG_level_set_swc_tools

#include "level_set_conforming_refine.h"   // LSSeg, LSSegDist, LSLoadSWC
#include "lib_disc/function_spaces/grid_function.h"

namespace ug {
namespace LevelSet {

template <typename TFunction>
void InterpolateSWCPhi(SmartPtr<TFunction> u, const char* cmp, const char* swcFile,
                       number radiusOffset, bool planar)
{
	static const int dim = TFunction::dim;
	std::vector<LSSeg> seg;
	if (!LSLoadSWC(swcFile, seg, planar))
		UG_THROW("InterpolateSWCPhi: cannot load SWC '" << swcFile << "'");

	const size_t fct = u->fct_id_by_name(cmp);
	typename TFunction::domain_type::position_accessor_type& aaPos =
			u->domain()->position_accessor();
	typedef typename TFunction::template traits<Vertex>::const_iterator vrt_iter;

	for (vrt_iter it = u->template begin<Vertex>(); it != u->template end<Vertex>(); ++it)
	{
		Vertex* v = *it;
		const MathVector<dim>& pp = aaPos[v];
		vector3 x(0, 0, 0);
		for (int d = 0; d < dim; ++d) x[d] = pp[d];
		number phi = 1e30;
		for (size_t k = 0; k < seg.size(); ++k) {
			const number d = LSSegDist(x, seg[k]) - radiusOffset;
			if (d < phi) phi = d;
		}
		std::vector<DoFIndex> ind;
		u->inner_dof_indices(v, fct, ind);
		if (!ind.empty()) DoFRef(*u, ind[0]) = phi;
	}
	u->set_storage_type(PST_CONSISTENT);   // values set per-vertex by hand -> mark consistent (VTK/markers)
}

//	Fill dst := sign_eps(src) = src / sqrt(src^2 + eps^2) (smoothed sign), per vertex. Used to build
//	the frozen reinitialisation source sign(phi0) as a GridFunction added to the RHS by operator
//	splitting (FV1_Convectionhang::set_source only takes a constant, not a field).
template <typename TFunction>
void FillSignEps(SmartPtr<TFunction> dst, SmartPtr<TFunction> src, const char* cmp, number eps)
{
	const size_t fct = dst->fct_id_by_name(cmp);
	typedef typename TFunction::template traits<Vertex>::const_iterator vrt_iter;
	for (vrt_iter it = dst->template begin<Vertex>(); it != dst->template end<Vertex>(); ++it)
	{
		Vertex* v = *it;
		std::vector<DoFIndex> is, id;
		src->inner_dof_indices(v, fct, is);
		dst->inner_dof_indices(v, fct, id);
		if (is.empty() || id.empty()) continue;
		const number p = DoFRef(*src, is[0]);
		DoFRef(*dst, id[0]) = p / std::sqrt(p * p + eps * eps);
	}
	dst->set_storage_type(PST_CONSISTENT);
}

//	Operator-split reinit source: u += dtau * sign_eps(phi0), per vertex. sign frozen at phi0 so the
//	interface (phi=0) is preserved. Replaces the missing Lua VecScaleAdd for the reinit pseudo-step.
template <typename TFunction>
void AddSignSource(SmartPtr<TFunction> u, SmartPtr<TFunction> phi0, const char* cmp,
                   number dtau, number eps)
{
	const size_t fct = u->fct_id_by_name(cmp);
	typedef typename TFunction::template traits<Vertex>::const_iterator vrt_iter;
	for (vrt_iter it = u->template begin<Vertex>(); it != u->template end<Vertex>(); ++it)
	{
		Vertex* v = *it;
		std::vector<DoFIndex> iu, i0;
		u->inner_dof_indices(v, fct, iu);
		phi0->inner_dof_indices(v, fct, i0);
		if (iu.empty() || i0.empty()) continue;
		const number p = DoFRef(*phi0, i0[0]);
		DoFRef(*u, iu[0]) += dtau * (p / std::sqrt(p * p + eps * eps));
	}
	u->set_storage_type(PST_CONSISTENT);
}

}  // namespace LevelSet
}  // namespace ug

#endif
