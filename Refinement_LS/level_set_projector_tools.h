/*
 * Templated bridge: read a level-set GridFunction φ into a LevelSetProjector's
 * per-vertex attachments (nodal φ + recovered nodal ∇φ). This is the lib_disc →
 * lib_grid bridge (the projector is non-templated lib_grid; the GridFunction is
 * templated lib_disc). Registered per domain-algebra.
 *
 * Workflow: build domain → interpolate/solve φ → SetLevelSetProjectorLSF(proj, u,
 * "phi") → attach proj to the refiner → refine.
 */
#ifndef __H__UG_level_set_projector_tools
#define __H__UG_level_set_projector_tools

#include "level_set_projector.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug {
namespace LevelSet {

template <typename TFunction>
void SetLevelSetProjectorLSF(SmartPtr<LevelSetProjector> proj,
                             SmartPtr<TFunction> u, const char* cmp)
{
	static const int dim = TFunction::dim;
	typedef typename TFunction::element_type element_type;   // Face(2D) / Volume(3D)
	typedef typename TFunction::template traits<Vertex>::const_iterator vrt_iter;
	typedef typename TFunction::template traits<element_type>::const_iterator elem_iter;

	MultiGrid& mg = *u->domain()->grid();
	proj->attach_grid(mg);

	const size_t fct = u->fct_id_by_name(cmp);
	typename TFunction::domain_type::position_accessor_type& aaPos =
			u->domain()->position_accessor();

//	temporary nodal-gradient accumulation (+ area weight)
	Attachment<vector3> aG;   Attachment<number> aW;
	mg.attach_to_vertices(aG); mg.attach_to_vertices(aW);
	Grid::VertexAttachmentAccessor<Attachment<vector3> > aaG(mg, aG);
	Grid::VertexAttachmentAccessor<Attachment<number> >  aaW(mg, aW);

//	1) nodal φ from the GridFunction; zero the gradient accumulators
	for (vrt_iter it = u->template begin<Vertex>(); it != u->template end<Vertex>(); ++it) {
		Vertex* v = *it;
		std::vector<DoFIndex> ind;
		u->inner_dof_indices(v, fct, ind);
		const number phiv = ind.empty() ? 0.0 : DoFRef(*u, ind[0]);
		proj->set_vertex_data(v, phiv, vector3(0, 0, 0));
		aaG[v] = vector3(0, 0, 0); aaW[v] = 0.0;
	}

//	2) recover nodal ∇φ = area-weighted average of the constant P1 element gradients
	MathMatrix<dim, dim> JTInv;
	std::vector<MathVector<dim> > vLocalGrad, vGlobalGrad;
	std::vector<MathVector<dim> > vCorner;

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
		const number elemSize = ElementSize<dim>(roid, &vCorner[0]);
		if (elemSize < 1e-14) continue;   // skip degenerate (e.g. constrained) elements
		map.update(&vCorner[0]);
		map.jacobian_transposed_inverse(JTInv, localIP);

	//	element-constant global gradient of φ
		MathVector<dim> gElem; VecSet(gElem, 0.0);
		for (size_t sh = 0; sh < numSH; ++sh) {
			MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);
			Vertex* vert = elem->vertex(sh);
			std::vector<DoFIndex> ind;
			u->inner_dof_indices(vert, fct, ind);
			VecScaleAppend(gElem, ind.empty() ? 0.0 : DoFRef(*u, ind[0]), vGlobalGrad[sh]);
		}
	//	accumulate to corners (weighted by element size)
		for (size_t sh = 0; sh < numSH; ++sh) {
			Vertex* vert = elem->vertex(sh);
			vector3 g3(0, 0, 0);
			for (int d = 0; d < dim; ++d) g3[d] = gElem[d];
			VecScaleAppend(aaG[vert], elemSize, g3);
			aaW[vert] += elemSize;
		}
	}

//	3) normalise and write the recovered ∇φ back onto the projector
	for (vrt_iter it = u->template begin<Vertex>(); it != u->template end<Vertex>(); ++it) {
		Vertex* v = *it;
		vector3 g = aaG[v];
		if (aaW[v] > 0.0) VecScale(g, g, 1.0 / aaW[v]);
		std::vector<DoFIndex> ind;
		u->inner_dof_indices(v, fct, ind);
		const number phiv = ind.empty() ? 0.0 : DoFRef(*u, ind[0]);
		proj->set_vertex_data(v, phiv, g);
	}

	mg.detach_from_vertices(aG);
	mg.detach_from_vertices(aW);
}

}  // namespace LevelSet
}  // namespace ug

#endif
