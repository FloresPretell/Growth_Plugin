/*
 * φ-driven ANISOTROPIC marker — the "orientation of refinement inside" half.
 *
 * The LevelSetProjector places new vertices (on φ=0); this marker chooses WHICH
 * edges of each band element split, guided by the level-set normal ∇φ, so the
 * refinement is anisotropic and oriented to the level set:
 *   acrossInterface=true  → split edges ‖ ∇φ  → thin cells ACROSS φ=0 (fine radial
 *                            band; what the interface/flux wants);
 *   acrossInterface=false → split edges ⊥ ∇φ  → fine ALONG the interface (tangential).
 *
 * Marks each band element RM_ANISOTROPIC(=RM_CLOSURE) + selected edges RM_REFINE,
 * then the caller runs refiner:refine(). Same mechanism as NeuriteRefMarkAdjuster,
 * driven by ∇φ instead of a neurite axis. (Produces hanging nodes unless globally
 * consistent — fine for transport terms; keep cut cells conforming for the flux.)
 */
#ifndef __H__UG_level_set_marker_tools
#define __H__UG_level_set_marker_tools

#include <vector>
#include <cmath>
#include <map>
#include "lib_grid/refinement/refiner_interface.h"
#include "lib_grid/multi_grid.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/domain_util.h"

namespace ug {
namespace LevelSet {

template <typename TFunction>
size_t MarkAnisotropicLevelSet(SmartPtr<IRefiner> refiner, SmartPtr<TFunction> u,
                               const char* cmp, number band, bool acrossInterface,
                               bool interior = false, int maxLvl = -1)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;

	MultiGrid& mg = *u->domain()->grid();
	typename TFunction::domain_type::position_accessor_type& aaPos =
			u->domain()->position_accessor();
	const size_t fct = u->fct_id_by_name(cmp);
	IRefiner& ref = *refiner;

	MathMatrix<dim, dim> JTInv;
	std::vector<MathVector<dim> > vLocalGrad, vGlobalGrad, vCorner;
	size_t nMarked = 0;

	for (elem_iter it = u->template begin<element_type>();
	     it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
	//	level cap: never refine an element already at maxLvl (bounds depth across frames
	//	so a cell that stays in the moving band does not refine unboundedly). -1 = no cap.
		if (maxLvl >= 0 && (int)mg.get_level(elem) >= maxLvl) continue;
		const ReferenceObjectID roid = elem->reference_object_id();
		const LocalShapeFunctionSet<dim>& lsfs =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));
		DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid);
		const MathVector<dim> localIP = ReferenceElementCenter<dim>(roid);
		const size_t numSH = lsfs.num_sh();
		vLocalGrad.resize(numSH); vGlobalGrad.resize(numSH);
		lsfs.grads(&vLocalGrad[0], localIP);
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		if (ElementSize<dim>(roid, &vCorner[0]) < 1e-14) continue;  // skip degenerate/constrained
		map.update(&vCorner[0]);
		map.jacobian_transposed_inverse(JTInv, localIP);

	//	element-centre φ (band test) and element-constant ∇φ (orientation)
		number phiC = 0;
		MathVector<dim> gElem; VecSet(gElem, 0.0);
		for (size_t sh = 0; sh < numSH; ++sh) {
			MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(elem->vertex(sh), fct, ind);
			const number val = ind.empty() ? 0.0 : DoFRef(*u, ind[0]);
			phiC += val / (number)numSH;
			VecScaleAppend(gElem, val, vGlobalGrad[sh]);
		}
	//	interior mode: refine the WHOLE interior (φ≤band) so radial shells reach the
	//	deep core, not just the symmetric |φ|≤band shell around φ=0. This one-sided
	//	gate is what makes the refinement follow the tube instead of leaving an
	//	unrefined core.
		if (interior) { if (phiC > band) continue; }
		else          { if (band > 0 && std::fabs(phiC) > band) continue; }

	//	ANISOTROPIC directional split (∇φ edge selection) is only valid on
	//	quad/hex elements — for triangles/tets UG4 has no clean anisotropic
	//	refinement (Part-5 finding), and the hanging-node refiner asserts. So
	//	only quads/hexes take the directional path; simplices fall back to a
	//	clean ISOTROPIC (full) refinement of the band element.
		const bool tensorElem = (numSH == 4 && dim == 2) || (numSH == 8 && dim == 3);
		if (!acrossInterface || !tensorElem) {
			ref.mark(elem, RM_REFINE);
			++nMarked;
			continue;
		}

		vector3 n(0, 0, 0);
		for (int d = 0; d < dim; ++d) n[d] = gElem[d];
		const number nlen = VecLength(n);
		if (nlen < 1e-14) { ref.mark(elem, RM_REFINE); ++nMarked; continue; }
		VecScale(n, n, 1.0 / nlen);

	//	quad/hex: anisotropic (closure) element + edges ‖ ∇φ split (thin across φ=0)
		ref.mark(elem, RM_ANISOTROPIC);
		++nMarked;
		typename MultiGrid::traits<Edge>::secure_container edges;
		mg.associated_elements(edges, elem);
		for (size_t i = 0; i < edges.size(); ++i) {
			Edge* e = edges[i];
			vector3 p0(0, 0, 0), p1(0, 0, 0);
			for (int d = 0; d < dim; ++d) { p0[d] = aaPos[e->vertex(0)][d]; p1[d] = aaPos[e->vertex(1)][d]; }
			vector3 ed; VecSubtract(ed, p0, p1);
			const number el = VecLength(ed);
			if (el < 1e-14) continue;
			VecScale(ed, ed, 1.0 / el);
			const number scp = std::fabs(VecDot(ed, n));       // |cos(edge, ∇φ)|
			ref.mark(e, scp > 0.707 ? RM_REFINE : RM_ANISOTROPIC);   // split ‖ ∇φ
		}
	//	3D: the FACES must also be marked RM_ANISOTROPIC so the hanging-node refiner
	//	keeps face closure consistent across LEVELS (without this the 2nd anisotropic
	//	pass asserts in refine_face_with_normal_vertex — face is neither adaptive nor
	//	full). The face's own edge marks (set above, shared) decide how it splits.
	//	Mirrors NeuriteRefMarkAdjuster::change_face_mark.
		if (dim == 3) {
			typename MultiGrid::traits<Face>::secure_container faces;
			mg.associated_elements(faces, elem);
			for (size_t i = 0; i < faces.size(); ++i)
				ref.mark(faces[i], RM_ANISOTROPIC);
		}
	}
	return nMarked;
}

//	SOLID-TUBE marker: refine every element whose centre φ ≤ outerBand — i.e. the
//	whole INTERIOR of the level set (φ<0) plus a thin outer band. This is the
//	"refine inside" case (a solid refined tube like the neurite), as opposed to
//	MarkAnisotropicLevelSet which refines only the |φ|-band shell. Isotropic
//	(RM_REFINE), robust on simplices; composes with the LevelSetProjector which
//	rounds the interface. Returns the number of marked elements.
template <typename TFunction>
size_t MarkLevelSetTube(SmartPtr<IRefiner> refiner, SmartPtr<TFunction> u,
                        const char* cmp, number outerBand)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;

	const size_t fct = u->fct_id_by_name(cmp);
	typename TFunction::domain_type::position_accessor_type& aaPos =
			u->domain()->position_accessor();
	IRefiner& ref = *refiner;
	std::vector<MathVector<dim> > vCorner;
	size_t nMarked = 0;

	for (elem_iter it = u->template begin<element_type>();
	     it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		if (ElementSize<dim>(roid, &vCorner[0]) < 1e-14) continue;  // degenerate/constrained

		const size_t n = elem->num_vertices();
		number phiC = 0;
		for (size_t i = 0; i < n; ++i) {
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(elem->vertex(i), fct, ind);
			phiC += (ind.empty() ? 0.0 : DoFRef(*u, ind[0])) / (number)n;
		}
		if (phiC <= outerBand) { ref.mark(elem, RM_REFINE); ++nMarked; }
	}
	return nMarked;
}

// ---------------------------------------------------------------------------
//  Growth / moving-interface marker. As the level set advances (the neurite
//  grows), NEW interior volume appears that must be refined so the solver stays
//  accurate at the advancing front. Given the OLD field uOld and the NEW field
//  uNew (same approximation space, same component), mark exactly the elements the
//  interface just swept into: phiOld_centre > band AND phiNew_centre <= band, plus
//  a collar of width `collar` on the new-interior side. Work is proportional to the
//  swept shell, not the whole domain -- the second projector Nicole asked for.
//  Refine, then re-attach the LevelSetProjector from uNew so the new vertices land
//  on the advanced interface.
// ---------------------------------------------------------------------------
template <typename TFunction>
size_t MarkLevelSetGrowthBand(SmartPtr<IRefiner> refiner,
                              SmartPtr<TFunction> uOld, SmartPtr<TFunction> uNew,
                              const char* cmp, number band, number collar)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;

	const size_t fctOld = uOld->fct_id_by_name(cmp);
	const size_t fctNew = uNew->fct_id_by_name(cmp);
	typename TFunction::domain_type::position_accessor_type& aaPos =
			uNew->domain()->position_accessor();
	IRefiner& ref = *refiner;
	std::vector<MathVector<dim> > vCorner;
	size_t nMarked = 0;

	for (elem_iter it = uNew->template begin<element_type>();
	     it != uNew->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		if (ElementSize<dim>(roid, &vCorner[0]) < 1e-14) continue;  // degenerate/constrained

		const size_t n = elem->num_vertices();
		number phiOldC = 0, phiNewC = 0;
		for (size_t i = 0; i < n; ++i) {
			std::vector<DoFIndex> iOld, iNew;
			uOld->inner_dof_indices(elem->vertex(i), fctOld, iOld);
			uNew->inner_dof_indices(elem->vertex(i), fctNew, iNew);
			phiOldC += (iOld.empty() ? 0.0 : DoFRef(*uOld, iOld[0])) / (number)n;
			phiNewC += (iNew.empty() ? 0.0 : DoFRef(*uNew, iNew[0])) / (number)n;
		}
	//	swept into the interior this step: was outside/at the old band, now inside
	//	the new band (+ a collar so the fresh front is resolved, not just the shell)
		if (phiOldC > band && phiNewC <= band + collar) {
			ref.mark(elem, RM_REFINE); ++nMarked;
		}
	}
	return nMarked;
}

// ---------------------------------------------------------------------------
//  COARSEN marker (C++, fast — replaces the slow per-element Lua CoarsenCB in the
//  growth loop): mark every element whose centre |φ| >= wc for RM_COARSEN. The
//  multigrid refiner restricts to complete surface families, so only previously-
//  refined cells actually coarsen. Reads φ from the GridFunction (no per-element SWC
//  loop) -> orders of magnitude faster than the Lua callback.
// ---------------------------------------------------------------------------
template <typename TFunction>
size_t MarkLevelSetCoarsen(SmartPtr<IRefiner> refiner, SmartPtr<TFunction> u,
                           const char* cmp, number wc)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;
	const size_t fct = u->fct_id_by_name(cmp);
	typename TFunction::domain_type::position_accessor_type& aaPos =
			u->domain()->position_accessor();
	IRefiner& ref = *refiner;
	std::vector<MathVector<dim> > vCorner;
	size_t nMarked = 0;
	for (elem_iter it = u->template begin<element_type>();
	     it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		if (ElementSize<dim>(roid, &vCorner[0]) < 1e-14) continue;
		const size_t n = elem->num_vertices();
		number phiC = 0;
		for (size_t i = 0; i < n; ++i) {
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(elem->vertex(i), fct, ind);
			phiC += (ind.empty() ? 0.0 : DoFRef(*u, ind[0])) / (number)n;
		}
		if (std::fabs(phiC) >= wc) { ref.mark(elem, RM_COARSEN); ++nMarked; }
	}
	return nMarked;
}

//	VERIFICATION: count surface elements CUT by the interface (min vertex φ < 0 < max vertex φ)
//	whose multigrid level is BELOW maxLvl — i.e. interface cells the adaptivity has NOT fully
//	refined. Returns 0 iff every interface cell sits at the finest level: proof that the
//	adaptive band keeps up with the moving front (no φ=0 crossing left in a coarse cell).
template <typename TFunction>
size_t CountInterfaceBelowLevel(SmartPtr<TFunction> u, const char* cmp, int maxLvl)
{
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;
	MultiGrid& mg = *u->domain()->grid();
	const size_t fct = u->fct_id_by_name(cmp);
	size_t nBelow = 0;
	for (elem_iter it = u->template begin<element_type>();
	     it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const size_t n = elem->num_vertices();
		number pmin = 1e30, pmax = -1e30;
		for (size_t i = 0; i < n; ++i) {
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(elem->vertex(i), fct, ind);
			const number v = ind.empty() ? 0.0 : DoFRef(*u, ind[0]);
			if (v < pmin) pmin = v;
			if (v > pmax) pmax = v;
		}
		if (pmin < 0.0 && pmax > 0.0 && (int)mg.get_level(elem) < maxLvl) ++nBelow;
	}
	return nBelow;
}

//	VERIFICATION for reinitialisation: max over band elements (|phi_centre|<=band) of ||∇φ|−1|.
//	0 == φ is a perfect signed distance in the band; large == φ drifted (needs reinit). Uses the
//	same element-constant ∇φ recovery as MarkAnisotropicLevelSet. band<=0 => whole domain.
template <typename TFunction>
number MaxGradDevInBand(SmartPtr<TFunction> u, const char* cmp, number band)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;
	typename TFunction::domain_type::position_accessor_type& aaPos = u->domain()->position_accessor();
	const size_t fct = u->fct_id_by_name(cmp);
	MathMatrix<dim, dim> JTInv;
	std::vector<MathVector<dim> > vLocalGrad, vGlobalGrad, vCorner;
	number maxDev = 0.0;
	for (elem_iter it = u->template begin<element_type>();
	     it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const ReferenceObjectID roid = elem->reference_object_id();
		const LocalShapeFunctionSet<dim>& lsfs =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));
		DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid);
		const MathVector<dim> localIP = ReferenceElementCenter<dim>(roid);
		const size_t numSH = lsfs.num_sh();
		vLocalGrad.resize(numSH); vGlobalGrad.resize(numSH);
		lsfs.grads(&vLocalGrad[0], localIP);
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		if (ElementSize<dim>(roid, &vCorner[0]) < 1e-14) continue;
		map.update(&vCorner[0]);
		map.jacobian_transposed_inverse(JTInv, localIP);
		number phiC = 0; MathVector<dim> g; VecSet(g, 0.0);
		for (size_t sh = 0; sh < numSH; ++sh) {
			MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(elem->vertex(sh), fct, ind);
			const number val = ind.empty() ? 0.0 : DoFRef(*u, ind[0]);
			phiC += val / (number)numSH;
			VecScaleAppend(g, val, vGlobalGrad[sh]);
		}
		if (band > 0 && std::fabs(phiC) > band) continue;
		const number dev = std::fabs(VecLength(g) - 1.0);
		if (dev > maxDev) maxDev = dev;
	}
	return maxDev;
}

//	MEAN ||∇φ|−1| over band elements (robust companion to MaxGradDevInBand: the max is dominated
//	by a single pathological cell — e.g. a junction/hanging cell — so the mean better reflects
//	whether reinit restored the signed-distance property across the band).
template <typename TFunction>
number MeanGradDevInBand(SmartPtr<TFunction> u, const char* cmp, number band)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;
	typename TFunction::domain_type::position_accessor_type& aaPos = u->domain()->position_accessor();
	const size_t fct = u->fct_id_by_name(cmp);
	MathMatrix<dim, dim> JTInv;
	std::vector<MathVector<dim> > vLocalGrad, vGlobalGrad, vCorner;
	number sumDev = 0.0; size_t n = 0;
	for (elem_iter it = u->template begin<element_type>();
	     it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const ReferenceObjectID roid = elem->reference_object_id();
		const LocalShapeFunctionSet<dim>& lsfs =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));
		DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid);
		const MathVector<dim> localIP = ReferenceElementCenter<dim>(roid);
		const size_t numSH = lsfs.num_sh();
		vLocalGrad.resize(numSH); vGlobalGrad.resize(numSH);
		lsfs.grads(&vLocalGrad[0], localIP);
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		if (ElementSize<dim>(roid, &vCorner[0]) < 1e-14) continue;
		map.update(&vCorner[0]);
		map.jacobian_transposed_inverse(JTInv, localIP);
		number phiC = 0; MathVector<dim> g; VecSet(g, 0.0);
		for (size_t sh = 0; sh < numSH; ++sh) {
			MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(elem->vertex(sh), fct, ind);
			const number val = ind.empty() ? 0.0 : DoFRef(*u, ind[0]);
			phiC += val / (number)numSH;
			VecScaleAppend(g, val, vGlobalGrad[sh]);
		}
		if (band > 0 && std::fabs(phiC) > band) continue;
		sumDev += std::fabs(VecLength(g) - 1.0); ++n;
	}
	return (n > 0) ? sumDev / (number)n : 0.0;
}

//	One EXPLICIT nodal reinitialisation step:  phi_i += dtau * sign_eps(phi0_i) * (1 - |∇phi|_i),
//	with |∇phi|_i the area-weighted nodal gradient magnitude. This is the correct pseudo-time
//	eikonal update (∂phi/∂tau = sign(phi0)(1-|∇phi|)); done NODALLY in C++ because a plain
//	convection-upwind (FV1_Convectionhang) is unstable for this Hamilton-Jacobi equation. The
//	frozen sign_eps(phi0) vanishes at the interface, so phi=0 is preserved (to O(eps)). Explicit:
//	use dtau < h_fine. Call n times for a sweep. band<=0 => whole domain, else only |phi0|<=band.
template <typename TFunction>
void ReinitStepNodal(SmartPtr<TFunction> u, SmartPtr<TFunction> phi0, const char* cmp,
                     number dtau, number eps, number band)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;
	typename TFunction::domain_type::position_accessor_type& aaPos = u->domain()->position_accessor();
	const size_t fct = u->fct_id_by_name(cmp);

	std::map<Vertex*, MathVector<dim> > gsum;   // area-weighted sum of element ∇phi
	std::map<Vertex*, number>           wsum;   // sum of weights (element volumes)
	MathMatrix<dim, dim> JTInv;
	std::vector<MathVector<dim> > vLocalGrad, vGlobalGrad, vCorner;

	//	pass 1: element-constant ∇phi -> accumulate onto vertices (reads the OLD u only)
	for (elem_iter it = u->template begin<element_type>(); it != u->template end<element_type>(); ++it)
	{
		element_type* elem = *it;
		const ReferenceObjectID roid = elem->reference_object_id();
		const LocalShapeFunctionSet<dim>& lsfs =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));
		DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid);
		const MathVector<dim> localIP = ReferenceElementCenter<dim>(roid);
		const size_t numSH = lsfs.num_sh();
		vLocalGrad.resize(numSH); vGlobalGrad.resize(numSH);
		lsfs.grads(&vLocalGrad[0], localIP);
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		const number vol = ElementSize<dim>(roid, &vCorner[0]);
		if (vol < 1e-14) continue;
		map.update(&vCorner[0]);
		map.jacobian_transposed_inverse(JTInv, localIP);
		MathVector<dim> g; VecSet(g, 0.0);
		for (size_t sh = 0; sh < numSH; ++sh) {
			MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);
			std::vector<DoFIndex> ind; u->inner_dof_indices(elem->vertex(sh), fct, ind);
			VecScaleAppend(g, ind.empty() ? 0.0 : DoFRef(*u, ind[0]), vGlobalGrad[sh]);
		}
		for (size_t sh = 0; sh < numSH; ++sh) {
			Vertex* v = elem->vertex(sh);
			VecScaleAppend(gsum[v], vol, g);
			wsum[v] += vol;
		}
	}

	//	pass 2: nodal update phi_i += dtau * sign_eps(phi0_i) * (1 - |∇phi|_i)
	typedef typename TFunction::template traits<Vertex>::const_iterator vrt_iter;
	for (vrt_iter it = u->template begin<Vertex>(); it != u->template end<Vertex>(); ++it)
	{
		Vertex* v = *it;
		typename std::map<Vertex*, number>::iterator wi = wsum.find(v);
		if (wi == wsum.end() || wi->second < 1e-30) continue;
		std::vector<DoFIndex> iu, i0;
		u->inner_dof_indices(v, fct, iu);
		phi0->inner_dof_indices(v, fct, i0);
		if (iu.empty() || i0.empty()) continue;
		const number p0 = DoFRef(*phi0, i0[0]);
		if (band > 0 && std::fabs(p0) > band) continue;
		MathVector<dim> gv; VecScale(gv, gsum[v], 1.0 / wi->second);
		const number gnorm = VecLength(gv);
		const number s = p0 / std::sqrt(p0 * p0 + eps * eps);
		DoFRef(*u, iu[0]) += dtau * s * (1.0 - gnorm);
	}
	u->set_storage_type(PST_CONSISTENT);
}

}  // namespace LevelSet
}  // namespace ug

#endif
