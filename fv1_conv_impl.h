/*
 * Author: Dmitry Logashenko
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 * Implementation of the FV-discretization of a convection equation (in a non-divergent form).
 */

/* UG4 headers: */
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

namespace ug{
namespace LevelSet{

////////////////////////////////////////////////////////////////////////////////
///	class constructor
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
FV1_Convection<TDomain>::
FV1_Convection
(
	SmartPtr<IConvectionShapes<dim> > upwind,
	const char * functions,
	const char * subsets
)
:	IElemDisc<TDomain> (functions, subsets),
	m_spUpwind (upwind), m_source (0), m_diffusion (0)
{
//	check number of functions
	if (this->num_fct () != 1)
		UG_THROW ("Wrong number of functions: The ElemDisc 'FV1_Convection'"
					" needs exactly 1 symbolic function.");
	
//	no derivatives of the velocity
	m_imVelocity.set_comp_lin_defect (false);
	
//	register imports
	this->register_import (m_imVelocity);

//	register assemble functions
	register_all_loc_discr_funcs ();
}
	
////////////////////////////////////////////////////////////////////////////////
//	check the grid and the shape functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void FV1_Convection<TDomain>::prepare_setting
(
	const std::vector<LFEID> & vLfeID,
	bool bNonRegular
)
{
//	check the grid
	if (bNonRegular)
		UG_THROW ("ERROR in FV1_Convection:"
				" The discretization does not support hanging nodes.\n");

//	check number of the components
	if (vLfeID.size () != 1)
		UG_THROW ("FV1_Convection: Only one component is supported.");

//	check that these are the Nedelec elements
	if (vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW ("FV1_Convection: This discretization works with the piecewise linear shape functions only.");
}

////////////////////////////////////////////////////////////////////////////////
// assembling
////////////////////////////////////////////////////////////////////////////////

/// prepares the loop over the elements: checks whether the parameters are set, ...
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::prepare_element_loop
(
	ReferenceObjectID roid, ///< only elements with this roid are looped over
	int si ///< and only in this subdomain
)
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static const int refDim = TElem::dim;
	
//	get the FV geometry
	static TFVGeom& geo = GeomProvider<TFVGeom>::get ();
	
//	check the velocity and the upwind method
	if(! m_imVelocity.data_given ())
		UG_THROW("FV1_Convection: No velocity assigned.");
	if(m_spUpwind.invalid ())
		UG_THROW("FV1_Convection: Upwind has not been set.");
	
//	set the local SCVF integration points
	const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	m_imVelocity.template set_local_ips<refDim> (vSCVFip, numSCVFip, false);
	
//	set the FV geometry type in the upwind method
	m_spUpwind->template set_geometry_type<TFVGeom> (geo);
}

/// finalizes the loop over the elements: clear the source
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::finish_element_loop ()
{
}

/// prepares a given element for assembling: computes the discretization of the rot-rot operator
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::prepare_element
(
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	ReferenceObjectID roid,  // id of reference element used for assembling
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
//	update the FV geometry
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try
	{
		geo.update (elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("FV1_Convection: Cannot update the Finite Volume Geometry.");
	
//	set the local SCVF integration points
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	m_imVelocity.set_global_ips (vSCVFip, numSCVFip);
}

/// assembles the local stiffness matrix
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::ass_JA_elem
(
	LocalMatrix & J,
	const LocalVector & u,
	GridObject * elem,
	const position_type vCornerCoords []
)
{
//	get the FV geometry
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static TFVGeom& geo = GeomProvider<TFVGeom>::get ();

//	update the upwind method
	if(! m_spUpwind->update (&geo, m_imVelocity.values (), NULL, false))
		UG_THROW("ERROR in 'FV1_Convection: Cannot compute convection shapes.\n");
	const size_t numConvShapes = m_spUpwind->num_sh();
	
// 	loop SCVFs
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf (ip);

		for(size_t sh = 0; sh < numConvShapes; ++sh)
		{
			const number D_conv_flux = (* m_spUpwind) (ip, sh);
			
		//	fluxes through the SCVF (due the convection term in the divergence form)
			J(_U_, scvf.from(), _U_, sh) += D_conv_flux;
			J(_U_, scvf.to(),   _U_, sh) -= D_conv_flux;
			
		//	sink due to the divergence
			J(_U_, scvf.from(), _U_, scvf.from()) -= D_conv_flux;
			J(_U_, scvf.to(),   _U_, scvf.to()  ) += D_conv_flux;
		}
	}
	
//	assemble the diffusion
	if(m_diffusion != 0)
	{
	//	loop SCVFs
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf (ip);
			
		// 	loop shape functions
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{
				const number D_diff_flux = m_diffusion * VecDot(scvf.global_grad(sh), scvf.normal());
				J(_U_, scvf.from(), _U_, sh) -= D_diff_flux;
				J(_U_, scvf.to()  , _U_, sh) += D_diff_flux;
			}
		}
	}
}

/// assembles the local defect
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::ass_dA_elem
(
	LocalVector & d,
	const LocalVector & u,
	GridObject * elem,
	const position_type vCornerCoords []
)
{
//	get the FV geometry
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	update the upwind method
	if(! m_spUpwind->update (&geo, m_imVelocity.values (), NULL, false))
		UG_THROW("ERROR in 'FV1_Convection: Cannot compute convection shapes.\n");
	const size_t numConvShapes = m_spUpwind->num_sh();
	
// 	loop SCVFs
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf (ip);

	//	fluxes through the SCVF (due the convection term in the divergence form)
		number conv_flux = 0;
		for(size_t sh = 0; sh < numConvShapes; ++sh)
			conv_flux += u (_U_, sh) * (* m_spUpwind) (ip, sh);
		d(_U_, scvf.from()) += conv_flux;
		d(_U_, scvf.to()  ) -= conv_flux;
			
	//	sink due to the divergence
		for(size_t sh = 0; sh < numConvShapes; ++sh)
		{
			d(_U_, scvf.from()) -= u (_U_, scvf.from()) * (* m_spUpwind) (ip, sh);
			d(_U_, scvf.to()  ) += u (_U_, scvf.to()  ) * (* m_spUpwind) (ip, sh);
		}
	}
	
//	assemble the diffusion
	if(m_diffusion != 0)
	{
	//	loop SCVFs
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf (ip);
			
		//	to compute D \nabla c
			MathVector<dim> grad;

		// 	compute gradient at ip
			VecSet(grad, 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad, u(_U_, sh), scvf.global_grad(sh));

		// 	compute flux
			const number diff_flux = m_diffusion * VecDot(grad, scvf.normal());

		// 	add to local defect
			d(_U_, scvf.from()) -= diff_flux;
			d(_U_, scvf.to()  ) += diff_flux;
		}
	}
}

/// computes the mass matrix of a time-dependent problem
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::ass_JM_elem
(
	LocalMatrix & J, 
	const LocalVector & u,
	GridObject * elem,
	const position_type vCornerCoords []
)
{
//	get the FV geometry
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop SCVs
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv (ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		J(_U_, co, _U_, co) += scv.volume ();
	}
}

/// computes the mass part of the defect of a time-dependent problem
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::ass_dM_elem
(
	LocalVector & d,
	const LocalVector & u,
	GridObject * elem,
	const position_type vCornerCoords []
)
{
//	get the FV geometry
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop SCVs
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv (ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		d(_U_, co) += u (_U_, co) * scv.volume ();
	}
}

/// assembles the right-hand side
template<typename TDomain>
template<typename TElem>
void FV1_Convection<TDomain>::ass_rhs_elem
(
	LocalVector & b,
	GridObject * elem,
	const position_type vCornerCoords []
)
{
//	get the FV geometry
	typedef FV1Geometry<TElem, dim> TFVGeom;
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	
// 	loop SCVs
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv (ip);

	// 	get associated node
		const int co = scv.node_id();

	//	add the source
		b(_U_, co) += m_source * scv.volume ();
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assembling functions
////////////////////////////////////////////////////////////////////////////////

/// registers the local assembler functions for all the elements and dimensions
template<typename TDomain>
void FV1_Convection<TDomain>::register_all_loc_discr_funcs ()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList> (RegisterLocalDiscr (this));
}

/// registers the local assembler functions for a given element
template<typename TDomain>
template<typename TElem> // the element to register for
void FV1_Convection<TDomain>::register_loc_discr_func ()
{
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	
	this->clear_add_fct(id);
	
	this->set_prep_elem_loop_fct(id, & this_type::template prepare_element_loop<TElem>);
	this->set_prep_elem_fct		(id, & this_type::template prepare_element<TElem>);
	this->set_fsh_elem_loop_fct	(id, & this_type::template finish_element_loop<TElem>);
	this->set_add_jac_A_elem_fct(id, & this_type::template ass_JA_elem<TElem>);
	this->set_add_jac_M_elem_fct(id, & this_type::template ass_JM_elem<TElem>);
	this->set_add_def_A_elem_fct(id, & this_type::template ass_dA_elem<TElem>);
	this->set_add_def_M_elem_fct(id, & this_type::template ass_dM_elem<TElem>);
	this->set_add_rhs_elem_fct	(id, & this_type::template ass_rhs_elem<TElem>);
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */
