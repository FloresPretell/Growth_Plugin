/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Authors: Christian Wehner, Dmitry Logashenko
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
 * High-resolution flux-based level set method
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET__HR_FB_LSM_DISCR_H__
#define __H__UG__PLUGINS__LEVEL_SET__HR_FB_LSM_DISCR_H__

#include <string>

// ug4 headers
#include "common/common.h"
#include "lib_algebra/operator/debug_writer.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace LevelSet{

/**
 * Class of an explicit discretization of the Level-Set-Equation
 *
 * This class implements the hight-resolution flux-based level-set method
 * for the discretization of the level-set equation of the form
 * \f{eqnarray*}{
 *   u_t + \nabla \cdot (u \, \mathbf{V}) = u \, \nabla \mathbf{V} + r,
 * \f}
 * where
 * <ul>
 * <li> \f$ u \f$			(unknown, scalar) the level-set function (solution)
 * <li> \f$ \mathbf{V} \f$	(given, vector) advection velocity
 * <li> \f$ r \f$			(given, scalar) source term
 * </ul>
 * The advection velocity consists of two parts:
 * \f{eqnarray*}{
 *  \mathbf{V} = \delta \frac{\nabla u}{\| \nabla u \|} + \gamma \mathb{V}_{\mathrm{user}},
 * \f}
 * where
 * <ul>
 * <li> \f$ \delta \f$		(given, scalar) a scaling factor for the normal velocity
 * <li> \f$ \mathbf{V} \f$	(given, vector) user-given velocity vector field
 * <li> \f$ \gamma \f$		(given, scalar) a scaling factor for the user-given velocity
 * </ul>
 *
 * For the details of the method, cf.
 * <ul>
 * <li> P. Frolkovic, K. Mikula, High-Resolution Flux-Based Level Set Method, SIAM. J. Sci. Comput., Vol. 29(2), pp. 579--597, DOI: 10.1137/050646561
 * <li> P. Frolkovic, Immersed Interface Method For a Level Set Formulation of Problems With Moving Boundaries, Proceedings of ALGORITHMY 2012, pp. 32--41
 * </ul>
 *
 * \tparam	TGridFunction	grid function type
 */
template<typename TGridFunction>
class HiResFluxBasedLSM
:	public DebugWritingObject<typename TGridFunction::algebra_type>
{
///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

///	world dimension
	static const int dim = domain_type::dim;

///	grid type
	typedef typename domain_type::grid_type grid_type;

/// type of the position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	
///	type of gradient attachment
	typedef Attachment<MathVector<dim> > ADimVector;

///	type of volume attachment accessor
	typedef typename Grid::VertexAttachmentAccessor<ANumber> t_aaVol;

///	type of gradient attachment accessor
	typedef typename Grid::VertexAttachmentAccessor<ADimVector> t_aaGrad;
	
///	type of the attachment accessor for the updates
	typedef typename Grid::VertexAttachmentAccessor<ANumber> t_aaUpd;

///	type of the attachment accessor for corners of intersected elements
	typedef typename Grid::VertexAttachmentAccessor<ABool> t_aaCoIE;

/// type of base grid object
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object ElemType;

/// edge iterator
	typedef typename TGridFunction::template traits<Edge>::const_iterator EdgeConstIterator;

/// vertex base iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexConstIterator;
	
///	grid element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

///	max. number of corners and subcontrol volume faces
	static const size_t maxNumCo = DimFV1Geometry<dim>::maxNumSCV;
	static const size_t maxNumIP = DimFV1Geometry<dim>::maxNumSCVF;
	
///	threshold for the level-set-function
	static number lsf_threshold () {return 1e-8;}
	
public:

///	Constructor
	HiResFluxBasedLSM()
	:	m_bVerbose (false), m_dt (0), m_nrOfSteps (1),
		m_time_control (false), m_maxCFL (0.95), m_minCFL (0.85),
		m_gamma (1), m_delta (0), m_divFree (false),
		m_firstOrder (false), m_antiderivSrc (false), m_elem_vel_vec (false),
		m_limiter (false),
		m_time (0), m_CFL (0)
	{
		set_source (0);
	}

///	Destructor
	virtual ~HiResFluxBasedLSM() {};
	
//	Controls

///	set the output level
	void set_verbose (bool b) { m_bVerbose = b; }

///	set the grid functions for the old and new solutions
	void set_solutions
	(
		SmartPtr<TGridFunction> uOld, ///< at the old time step
		SmartPtr<TGridFunction> uNew ///< at the new time step
	)
	{
		m_oldSol = uOld; m_newSol = uNew;
	}

///	set level-set function to specify the interface
	void set_LSF (SmartPtr<TGridFunction> spLSF) { m_spLSF = spLSF; }
	
///	set the signed-distance function for the computation of the effective time step length at the interface
	void set_SDF (SmartPtr<TGridFunction> spSDF) { m_spSDF = spSDF; }
	
///	set the potential for the computation of the velocity
	void set_vel_potential (SmartPtr<TGridFunction> spVelPot) { m_spVelPot = spVelPot; }

///	set time step
	void set_dt (number dt) { m_dt = dt; }
	
///	get time step
	number get_dt () { return m_dt; }
	
///	switch the time control on
	void set_time_control (number minCFL, number maxCFL)
		{ m_time_control = true; m_minCFL = minCFL; m_maxCFL = maxCFL; }
		
/// switch the time control off
	void set_time_control_off () { m_time_control = false; }

/// set nr of time steps to perform
	void set_nr_of_steps (size_t n) { m_nrOfSteps = n; }
	
///	set the time argument
	void set_time (number t) { m_time = t; }
/// get the current time argument
	number get_time () { return m_time; }
///	get the max. CFL
	number get_CFL () { return m_CFL; }
	
///	set scaling factor for the user-given velocity (if it is given by the user data)
	void set_gamma (number gamma) { m_gamma = gamma; }
	
///	set scaling factor for the gradient in the velocity (if it is used for the velocity)
	void set_delta (number delta) { m_delta = delta; }
	
///	sets the divergence free flag
	void set_divfree (bool b) { m_divFree = b; }
	
///	switches the first order upwind method on/off
	void set_first_order (bool b) { m_firstOrder = b; }
	
///	sets the use of the antiderivative in the discretization of the source term
	void set_antideriv_src (bool b) { m_antiderivSrc = b; }
	
///	sets the interpretation of the vector velocity
	void set_elem_vel_vec (bool b) { m_elem_vel_vec = b; }
	
///	sets whether to use the slope limiter
	void set_limiter (bool b) { m_limiter = b; }
	
///	set a constant source term (both the values)
	void set_source (number val) { m_source_pos = val; m_source_neg = - val; }
	
///	set the source only for the subdomain with the negative values of the LSF
	void set_source_neg (number val) { m_source_neg = val; }

///	set the source only for the subdomain with the positive values of the LSF
	void set_source_pos (number val) { m_source_pos = val; }
	
///	set velocity vector field
	void set_velocity (SmartPtr<CplUserData<MathVector<dim>, dim> > vel) {m_imVelocity = vel;}
#ifdef UG_FOR_LUA
///	set velocity vector field as lua function
	void set_velocity (const char* fctName) { set_velocity (LuaUserDataFactory<MathVector<dim>,dim>::create (fctName)); }
#endif
	
///	set normal velocity field
	void set_normal_velocity (SmartPtr<CplUserData<number, dim> > vel) {m_imNormalVel = vel;}
#ifdef UG_FOR_LUA
///	set normal velocity field as lua function
	void set_normal_velocity (const char* fctName) { set_normal_velocity (LuaUserDataFactory<number,dim>::create (fctName)); }
#endif
	
///	sets the Dirichlet values at the interface
	void set_interface_data (SmartPtr<CplUserData<number,dim> > d) { m_imInterfaceVal = d; }
#ifdef UG_FOR_LUA
///	sets the Dirichlet values at the interface as lua function
	void set_interface_data (const char* fctName) { set_interface_data (LuaUserDataFactory<number,dim>::create (fctName)); }
#endif

///	set the Dirichlet values as a user data object
	void set_dirichlet_data (SmartPtr<CplUserData<number,dim> > d) { m_imDirichlet = d; }
///	set constant Dirichlet values
	void set_dirichlet_data (number val) { set_dirichlet_data (make_sp (new ConstUserNumber<dim> (val))); }
#ifdef UG_FOR_LUA
///	set the Dirichlet values as a lua function
	void set_dirichlet_data (const char* fctName) { set_dirichlet_data (LuaUserDataFactory<number,dim>::create (fctName)); }
#endif

/// boundary condition subset handling: dirichlet boundary
	void set_dirichlet_boundary (const char* subsets);

/// boundary condition subset handling: outflow boundary
	void set_outflow_boundary (const char* subsets);
	
///	prerapre the object for the computation of the SDF
	void prepare_for_SDF
	(
		SmartPtr<TGridFunction> uOld, ///< solution at the old time step
		SmartPtr<TGridFunction> uNew ///< solution at the new time step
	)
	{
		set_solutions (uOld, uNew);
		set_SDF (uOld);
		set_vel_potential (uOld);
		set_delta (1);
		set_gamma (0);
		set_source (1);
	}
	
//	Computation

///	computes the time steps of the discretization of the level-set equation
	void advect ();
	
///	compute the normal velocity using a user-data object
	void compute_normal_vel
	(
		SmartPtr<CplUserData<MathVector<dim>, dim> > spVelField,
		SmartPtr<TGridFunction> spNormVel
	);
		
///	appends a vertical vector the the normal velocity (using a user-data object)
	void append_vertical_to_normal_vel
	(
		SmartPtr<CplUserData<number, dim> > spNVelField,
		SmartPtr<TGridFunction> spNormVel
	);
	
///	appends a vertical vector the the normal velocity (using a user-data object)
	void append_vertical_to_normal_vel
	(
		SmartPtr<TGridFunction> spNVelGF,
		const char * cmp,
		SmartPtr<TGridFunction> spNormVel
	)
	{
		append_vertical_to_normal_vel
		(
			SmartPtr<CplUserData<number, dim> > (new GridFunctionNumberData<TGridFunction> (spNVelGF, cmp)),
			spNormVel
		);
	}
	
#ifdef UG_FOR_LUA
///	compute the normal velocity using a user-data object
	void compute_normal_vel
	(
		const char* fctName,
		SmartPtr<TGridFunction> spNormVel
	)
	{
		compute_normal_vel (LuaUserDataFactory<MathVector<dim>,dim>::create (fctName), spNormVel);
	}
#endif
	
///	compare the solution with a given level-set function
	void compare_lsf_with
	(
		SmartPtr<TGridFunction> spLSF2, ///< the second level-set function
		number eps ///< tolerance
	);
	
///	specifies a grid function to save the Courant number to
	void save_CourantNumber_to
	(
		SmartPtr<TGridFunction> spCN ///< a grid function for the Current numbers
	)
	{
		m_spCourant = spCN;
	}

private:

//	Auxiliary tools

/// slope limiter
	void limit_grad(TGridFunction& uOld, t_aaGrad& aaGradient);

///	computes the scvf-update of the solution in an element
	inline void sol_update
	(
		bool redOrder,
		const MathVector<dim>& ip,
		const MathVector<dim>& x_up,
		number u_up,
		const MathVector<dim>& grad_up,
		const MathVector<dim>& vel_up,
		const MathVector<dim>& x_down,
		number u_down,
		const MathVector<dim>& grad_down,
		const MathVector<dim>& vel_down,
		number& corr_up,
		number& curr_down,
		number& src_up,
		number& src_down
	);
///	computes the bf-update of the solution in an element
	inline void bnd_sol_update
	(
		bool redOrder,
		const MathVector<dim>& bip,
		const MathVector<dim>& x,
		number u,
		const MathVector<dim>& grad,
		const MathVector<dim>& vel,
		number& curr
	);
///	gets corner velocity and source
	inline void get_nodal_vel
	(
		ElemType* elem,
		MathVector<dim> coCoord[],
		DimFV1Geometry<dim>& geo,
		LocalVector& u,
		MathVector<dim> grad[],
		MathVector<dim> co_vel[],
		int lsf_sign
	);
///	assemble local contributions of one element
	void assemble_element
	(
		ElemType* elem,
		DimFV1Geometry<dim>& geo,
		domain_type& grid,
		LocalVector& uOld,
		t_aaGrad& aaGradient,
		t_aaGrad& aaVelGrad,
		t_aaVol& aaVolume,
		int sign,
		t_aaUpd& aaUpdate,
		t_aaUpd* aaSrc
	);
///	get the velocity for a given SCVF in an element intersected by the interface
	void get_scvf_vel_on_if
	(
		DimFV1Geometry<dim>& geo,
		const typename DimFV1Geometry<dim>::SCVF& scvf,
		number u[],
		MathVector<dim> grad[],
		number lsf[],
		MathVector<dim>& from_co_vel,
		number& from_flux,
		MathVector<dim>& to_co_vel,
		number& to_flux
	);
/// get the velocity for a given BF in an element intersected by the interface
	void get_bf_vel_on_if
	(
		DimFV1Geometry<dim>& geo,
		const typename DimFV1Geometry<dim>::BF& bf,
		number u[],
		MathVector<dim> grad[],
		number lsf[],
		MathVector<dim>& co_vel,
		number& flux
	);
///	assemble an element intersected by the interface
	int assemble_cut_element
	(
		ElemType* elem,
		DimFV1Geometry<dim>& geo,
		domain_type& domain,
		LocalVector& uOld,
		LocalVector& locLSF,
		LocalVector& locVelPot,
		t_aaGrad& aaGradient,
		t_aaGrad& aaVelGrad,
		t_aaVol& aaVolume,
		t_aaUpd& aaUpdate,
		t_aaUpd* aaSrc,
		CplUserData<number,dim> * if_val_data,
		int si
	);

/// compute CV volumes
	void compute_volumes
	(
		TGridFunction& u,
		DimFV1Geometry<dim>& geo,
		ANumber& aVolume,
		t_aaVol& aaVolume
	);
	
///	compute gradients in an element
	void compute_elem_grad
	(
		DimFV1Geometry<dim> & geo,
		number uValue [],
		MathVector<dim> co_grad [],
		number * lsf,
		CplUserData<number,dim> * if_val_data,
		int si
	);
/// compute gradients and volumes
	void compute_vertex_grad
	(
		TGridFunction& u,
		DimFV1Geometry<dim>& geo,
		t_aaVol& aaVolume,
		ADimVector& aGradient,
		t_aaGrad& aaGradient,
		TGridFunction * pLSF = NULL,
		CplUserData<number,dim> * if_val_data = NULL
	);
	
///	sign of the LSF
	inline int lsf_sign
	(
		size_t noc,
		number lsf []
	);
	
///	extrapolation by the LSF
	inline void extrapolate_by_lsf
	(
		const CplUserData<number,dim> * if_val_data,
		int si,
		DimFV1Geometry<dim> & geo,
		number sol[],
		number lsf[],
		size_t base,
		number ext[]
	);
	
///	mark corners at the interface
	void mark_CoIE
	(
		grid_type& grid,
		ABool& aCoIE,
		t_aaCoIE& aaCoIE
	);

///	assign Dirichlet values
	void assign_dirichlet
	(
		TGridFunction& numsol
	);

private:

//	Grid functions:
	
	SmartPtr<TGridFunction> m_oldSol; ///< solution at the old time step
	SmartPtr<TGridFunction> m_newSol; ///< computed solution at the new time step
	
	SmartPtr<TGridFunction> m_spLSF; ///< Level-Set Function data (if any)
	SmartPtr<TGridFunction> m_spSDF; ///< Signed-Distance Function (for computations of the eff. dt at the interface)
	SmartPtr<TGridFunction> m_spVelPot; ///< vert.-centred potential for the computation of the velocity (or SPNULL)
	
	SmartPtr<TGridFunction> m_spCourant; ///< a grid function for the local Courant numbers (optional)
	
//	Parameters of the method:

	bool m_bVerbose; ///< whether to print more details

	number m_dt; ///< current time step
	size_t m_nrOfSteps; ///< number of time steps to compute
	bool m_time_control; ///< whether to compute the appropriate time step
	number m_maxCFL; ///< max. allowed Courant number in a time step
	number m_minCFL; ///< min. allowed Courant number in a time step
	
	number m_gamma; ///< scaling factor for the user-given velocity (if it is given by the user data)
	number m_delta; ///< scaling factor for the gradient of the potential in the velocity (if it is used for the velocity)
	bool m_divFree; ///< if the velocity field is divergence free
	
	bool m_firstOrder; ///< if to use the classic (first order) upwind method
	bool m_antiderivSrc; ///< if to use the antiderivative to discretize the source
	
	bool m_elem_vel_vec; ///< for the normal velocity, take the elem.-centered velocity
	
	bool m_limiter; ///< whether to use the slope limiter
	
	SubsetGroup m_neumann_sg; ///< subsets with the Neumann BC
	SubsetGroup m_dirichlet_sg; ///< subsets with the Dirichlet BC
	
	SmartPtr<CplUserData<MathVector<dim>, dim> > m_imVelocity; ///< data import for the velocity field (if used)
	SmartPtr<CplUserData<number, dim> > m_imNormalVel; ///< data import for the normal velocity field (if used)
	SmartPtr<CplUserData<number,dim> > m_imDirichlet; ///< data import for the Dirichlet values
	SmartPtr<CplUserData<number,dim> > m_imInterfaceVal; ///< Dirichlet values at the interface
	
///	Values of the source for the positive and negative values of the LSF. (If no LSF, only the former is used.)
	number m_source_pos, m_source_neg;
	
//	Temporary and computed data
	
	number m_time; ///< current time
	number m_CFL; ///< max. Courant number achieved in all the computed steps
};

} // end namespace LevelSet
} // end namespace ug

// include implementation
#include "hrfblsm_discr_impl.h"

#endif /* __H__UG__PLUGINS__LEVEL_SET__HR_FB_LSM_DISCR_H__ */

/* End of File */
