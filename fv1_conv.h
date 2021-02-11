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
 * FV-discretization of a convection equation (in a non-divergent form).
 */

#ifndef __H__UG__PLUGINS__LEVEL_SET__FV1_CONV_H__
#define __H__UG__PLUGINS__LEVEL_SET__FV1_CONV_H__

// basic ug4 headers
#include "common/common.h"

// library-specific headers
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace LevelSet{

/**
 * The class implements the 1-st order FV discretization of the convection equation
 * (in the non-divergence form) with a diffusion term:
 * \f{eqnarray*}{
 *   u_t + \mathbf{V} \cdot \nabla u - \nabla \cdot (d \nabla u) = r,
 * \f}
 * where \f$\mathbf{V}\f$ is the convection velocity.
 * 
 * For the discretization, the equation is reformulated into the divergence form with
 * a sink term:
 * \f{eqnarray*}{
 *   u_t + \nabla \cdot (u \, \mathbf{V}) = u \, \nabla \mathbf{V} + r,
 * \f}
 *
 * \tparam	TDomain		Domain type
 */
template <typename TDomain>
class FV1_Convection
	: public IElemDisc<TDomain>
{
private:
///	abbreviation for the local solution
	static const size_t _U_ = 0;
	
///	base class type
	typedef IElemDisc<TDomain> base_type;

///	own type
	typedef FV1_Convection<TDomain> this_type;

///	domain type
	typedef typename base_type::domain_type domain_type;

///	position type
	typedef typename base_type::position_type position_type;
	
///	world dimension
	static const int dim = base_type::dim;
	
///	upwind method type
	typedef IConvectionShapes<dim> conv_shape_type;
	
public:

//---- Constructor ----

/// constructor
	FV1_Convection
	(
		SmartPtr<IConvectionShapes<dim> > upwind, ///< upwind method to use
		const char * functions, ///< symbolic function name
		const char * subsets ///< full-dim. subsets where to assemble
	);

//---- Settings: ----

///	sets the upwind for the transport equation
	void set_upwind (SmartPtr<IConvectionShapes<dim> > upwind)
	{
		m_spUpwind = upwind;
	}

///	sets the convection velocity
	void set_velocity (SmartPtr<CplUserData<MathVector<dim>, dim> > vel)
	{
		m_imVelocity.set_data (vel);
	}

#ifdef UG_FOR_LUA

///	sets the convection velocity from a LUA function specified by name
	void set_velocity(const char* fctName)
	{
		set_velocity(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
	}
	
///	sets the convection velocity from a LUA function
	void set_velocity(LuaFunctionHandle fct)
	{
		set_velocity(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
	}
	
#endif

///	sets the value of the source
	void set_source (number r)
	{
		m_source = r;
	}
	
///	sets the diffusion coefficient
	void set_diffusion (number d)
	{
		m_diffusion = d;
	}

//---- Local discretization interface: ----
private:
	
///	check type of the grid and the trial space
	virtual void prepare_setting
	(
		const std::vector<LFEID> & vLfeID,
		bool bNonRegular
	);

//---- Assembling functions: ----
	
	template <typename TElem>
	void prepare_element_loop(ReferenceObjectID roid, int si);

	template <typename TElem>
	void prepare_element(const LocalVector& u, GridObject* elem, ReferenceObjectID roid, const position_type vCornerCoords[]);

	template <typename TElem>
	void finish_element_loop();

	template <typename TElem>
	void ass_JA_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dA_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_JM_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dM_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_rhs_elem(LocalVector& d, GridObject* elem, const position_type vCornerCoords[]);

//---- Registration of the template functions: ----
private:
	
	void register_all_loc_discr_funcs();

	struct RegisterLocalDiscr
	{
			RegisterLocalDiscr(this_type * pThis) : m_pThis(pThis) {}
			this_type * m_pThis;
			template< typename TElem > void operator() (TElem &)
				{m_pThis->register_loc_discr_func<TElem> ();}
	};

	template <typename TElem>
	void register_loc_discr_func ();

private:

	SmartPtr<conv_shape_type> m_spUpwind; ///< upwind method
	
	DataImport<MathVector<dim>, dim > m_imVelocity; ///< data import for the velocity field (if used)
	
	number m_source; ///< the source term \f$r\f$
	
	number m_diffusion; ///< the diffusion coefficient
	
}; // end class 

} // end namespace LevelSet
} // end namespace ug

#include "fv1_conv_impl.h"

#endif // __H__UG__PLUGINS__LEVEL_SET__FV1_CONV_H__

/* End of File */
