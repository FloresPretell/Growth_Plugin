/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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


// ug4 headers
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"

// plugin headers
#include "level_set.h"
#include "level_set_user_data.h"
#include "ls_analytic.h"
//#include "ls_curvature2d.h"
#include "hrfblsm_discr.h"
#include "ls_volume.h"
#include "ls_init.h"
#include "fv1_conv.h"
#include "normvel_util.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace LevelSet{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
	
//	1st order upwinded discretization of the convection equation
	{
		typedef FV1_Convection<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("FV1_Convection").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(SmartPtr<IConvectionShapes<dim> >,const char*,const char*)>("Upwind#Function(s)#Subset(s)")
			.add_method("set_upwind", static_cast<void (T::*)(SmartPtr<IConvectionShapes<dim> >)>(&T::set_upwind), "Upwind method", "Upwind (no, full, ...)")
			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_velocity), "", "Velocity Field")
#endif
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "value")
			.add_method("set_diffusion", static_cast<void (T::*)(number)>(&T::set_diffusion), "", "value")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1_Convection", tag);
	}
}

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	// typedef typename TAlgebra::vector_type vector_type;
	// typedef typename TAlgebra::matrix_type matrix_type;
	// typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;
	static const int dim = TDomain::dim;

// 	FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		// typedef typename function_type::domain_type domain_type;
		string name = string("FV1LevelSetDisc").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_dt", &T::set_dt)
			.add_method("set_time", &T::set_time)
			.add_method("set_reinit", &T::set_reinit)
			.add_method("set_divfree", &T::set_divfree)
			.add_method("set_info", &T::set_info)
			.add_method("set_timestep_nr",&T::set_timestep_nr)
			.add_method("set_nr_of_steps", &T::set_nr_of_steps)
			.add_method("advect_lsf", &T::advect_lsf)
			.add_method("set_delta", &T::set_delta)
			.add_method("set_gamma", &T::set_gamma)
			.add_method("set_limiter", &T::set_limiter)
			.add_method("set_dirichlet_boundary",&T::set_dirichlet_boundary)
			.add_method("set_outflow_boundary",&T::set_outflow_boundary)
			.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif
			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(number)>(&T::set_velocity), "", "Vel_x")
			.add_method("set_velocity", static_cast<void (T::*)(number,number)>(&T::set_velocity), "", "Vel_x, Vel_y")
			.add_method("set_velocity", static_cast<void (T::*)(number,number,number)>(&T::set_velocity), "", "Vel_x, Vel_y, Vel_z")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
#endif
			.add_method("set_dirichlet_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_dirichlet_data), "", "Source")
			.add_method("set_dirichlet_data", static_cast<void (T::*)(number)>(&T::set_dirichlet_data), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_dirichlet_data", static_cast<void (T::*)(const char*)>(&T::set_dirichlet_data), "", "Source")
#endif
			.add_method("compute_normal", &T::compute_normal)
			.add_method("compute_dnormal", &T::compute_dnormal)
			.add_method("compute_ddnormal", &T::compute_ddnormal)
			.add_method("get_time", &T::get_time)
			.add_method("runtimetest", &T::runtimetest)
			.add_method("init_ls_subsets", &T::init_ls_subsets)
			.add_method("create_ls_subsets", &T::create_ls_subsets)
			.add_method("update_ls_subsets" ,&T::update_ls_subsets)
			.add_method("set_elements_active", static_cast<void (T::*)(int)>(&T::set_elements_active))
			.add_method("set_elements_active", static_cast<void (T::*)(int,int)>(&T::set_elements_active))
			.add_method("set_elements_active", static_cast<void (T::*)(int,int,int)>(&T::set_elements_active))
			.add_method("set_elements_inactive", static_cast<void (T::*)(int)>(&T::set_elements_inactive))
			.add_method("set_elements_inactive", static_cast<void (T::*)(int,int)>(&T::set_elements_inactive))
			.add_method("set_elements_inactive", static_cast<void (T::*)(int,int,int)>(&T::set_elements_inactive))
			.add_method("set_nodes_active", static_cast<void (T::*)(int)>(&T::set_nodes_active))
			.add_method("set_nodes_active", static_cast<void (T::*)(int,int)>(&T::set_nodes_active))
			.add_method("set_nodes_active", static_cast<void (T::*)(int,int,int)>(&T::set_nodes_active))
			.add_method("set_nodes_inactive", static_cast<void (T::*)(int)>(&T::set_nodes_inactive))
			.add_method("set_nodes_inactive", static_cast<void (T::*)(int,int)>(&T::set_nodes_inactive))
			.add_method("set_nodes_inactive", static_cast<void (T::*)(int,int,int)>(&T::set_nodes_inactive))
			.add_method("overwrite",static_cast<bool (T::*)(function_type&,function_type&,function_type&,int)>(&T::overwrite))
			.add_method("overwrite",static_cast<bool (T::*)(function_type&,number,function_type&,int)>(&T::overwrite))
			.add_method("compute_error", &T::compute_error)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "FV1LevelSetDisc", tag);
	}
	
	// level set analytic functions
	{
		typedef LevelSetAnalytic<function_type> T;
		// typedef typename function_type::domain_type domain_type;
		string name = string("LevelSetAnalytic").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_time", &T::set_time)
			.add_method("get_time", &T::get_time)
			.add_method("fill_v_vec",&T::fill_v_vec)
			.add_method("init_function", &T::init_function)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LevelSetAnalytic", tag);
	}
	
	// computation of the curvature in 2d
//	{
//		typedef LevelSetCurvature<function_type> T;
//		typedef typename function_type::domain_type domain_type;
//		string name = string("LevelSetCurvature").append(suffix);
//		reg.add_class_<T>(name, grp)
//			.add_constructor()
//			.add_method("compute_curvature", &T::computeElementCurvatureOnGrid2d)
//			.add_method("compute_curvature_from_sides", &T::computeElementCurvatureFromSides)
//			.add_method("exact_curvature", static_cast<void (T::*)(number)>(&T::exact_curvature), "", "Exact curvature")
//			.set_construct_as_smart_pointer(true);
//		reg.add_class_to_group(name, "LevelSetCurvature", tag);
//	}

	// level set user data
	{
		string name = string("LevelSetUserData").append(suffix);
		typedef LevelSetUserData<function_type> T;
		typedef CplUserData<number, dim> TBase;
		typedef INewtonUpdate TBase2;
		reg.add_class_<T, TBase,TBase2>(name, grp)
			.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<function_type>)>("Approximation space, grid function")
			.add_method("set_inside", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_inside_data), "", "Inside data")
			.add_method("set_inside", static_cast<void (T::*)(number)>(&T::set_inside_data), "", "Inside data")
		#ifdef UG_FOR_LUA
			.add_method("set_inside", static_cast<void (T::*)(const char*)>(&T::set_inside_data), "", "Inside data")
		#endif
			.add_method("set_outside", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_outside_data), "", "Outside data")
			.add_method("set_outside", static_cast<void (T::*)(number)>(&T::set_outside_data), "", "Outside data")
		#ifdef UG_FOR_LUA
			.add_method("set_outside", static_cast<void (T::*)(const char*)>(&T::set_outside_data), "", "Outside data")
		#endif
			.add_method("set_eval_type", static_cast<void (T::*)(int)>(&T::set_eval_type), "", "Evaluation type")
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LevelSetUserData", tag);
	}

	// level set user vector data
	{
			string name = string("LevelSetUserVectorData").append(suffix);
			typedef LevelSetUserVectorData<function_type> T;
			typedef CplUserData<MathVector<dim>, dim> TBase;
			typedef INewtonUpdate TBase2;
			reg.add_class_<T, TBase,TBase2>(name, grp)
				.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<function_type>)>("Approximation space, grid function")
				.add_method("set_inside", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_inside_data), "", "Source")
				.add_method("set_inside", static_cast<void (T::*)(number)>(&T::set_inside_data), "", "F_x")
				.add_method("set_inside", static_cast<void (T::*)(number,number)>(&T::set_inside_data), "", "F_x, F_y")
				.add_method("set_inside", static_cast<void (T::*)(number,number,number)>(&T::set_inside_data), "", "F_x, F_y, F_z")
			#ifdef UG_FOR_LUA
				.add_method("set_inside", static_cast<void (T::*)(const char*)>(&T::set_inside_data), "", "Source Vector")
			#endif
				.add_method("set_outside", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_outside_data), "", "Source")
				.add_method("set_outside", static_cast<void (T::*)(number)>(&T::set_outside_data), "", "F_x")
				.add_method("set_outside", static_cast<void (T::*)(number,number)>(&T::set_outside_data), "", "F_x, F_y")
				.add_method("set_outside", static_cast<void (T::*)(number,number,number)>(&T::set_outside_data), "", "F_x, F_y, F_z")
			#ifdef UG_FOR_LUA
				.add_method("set_outside", static_cast<void (T::*)(const char*)>(&T::set_outside_data), "", "Source Vector")
			#endif
				.add_method("set_eval_type", static_cast<void (T::*)(int)>(&T::set_eval_type), "", "Evaluation type")
			.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "LevelSetUserVectorData", tag);
		}

	// two phase flow source data
		{
			string name = string("CRTwoPhaseSource").append(suffix);
			typedef CRTwoPhaseSource<function_type> T;
			typedef CplUserData<MathVector<dim>, dim> TBase;
			typedef INewtonUpdate TBase2;
			reg.add_class_<T, TBase,TBase2>(name, grp)
				.template add_constructor<void (*)(SmartPtr<ApproximationSpace<TDomain> >,SmartPtr<function_type>)>("Approximation space, grid function")
				.add_method("set_density", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_density), "", "Set density")
				.add_method("set_density", static_cast<void (T::*)(number)>(&T::set_density), "", "Set density")
			#ifdef UG_FOR_LUA
				.add_method("set_density", static_cast<void (T::*)(const char*)>(&T::set_density), "", "Set density")
			#endif
				.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_source), "", "Source")
				.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "F_x")
				.add_method("set_source", static_cast<void (T::*)(number,number)>(&T::set_source), "", "F_x, F_y")
				.add_method("set_source", static_cast<void (T::*)(number,number,number)>(&T::set_source), "", "F_x, F_y, F_z")
			#ifdef UG_FOR_LUA
				.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source Vector")
			#endif
				.add_method("set_sigma", static_cast<void (T::*)(number)>(&T::set_sigma), "", "Set sigma")
				.add_method("set_gravitation", static_cast<void (T::*)(number)>(&T::set_gravitation), "", "Set gravitation constant")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "CRTwoPhaseSource", tag);
		}

// 	HiResFluxBasedLSM
	{
		typedef HiResFluxBasedLSM<function_type> T;
		typedef DebugWritingObject<TAlgebra> TBase;
		
		// typedef typename function_type::domain_type domain_type;
		string name = string("HiResFluxBasedLSM").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("advect", static_cast<void (T::*) ()>(&T::advect), "", "Compute the time steps of the method")
			.add_method("set_solutions", static_cast<void (T::*)(SmartPtr<function_type>,SmartPtr<function_type>)>(&T::set_solutions), "old#new", "Solutions at the time levels")
			.add_method("set_LSF", static_cast<void (T::*)(SmartPtr<function_type>)>(&T::set_LSF), "", "Level-Set function to specity the interface")
			.add_method("set_SDF", static_cast<void (T::*)(SmartPtr<function_type>)>(&T::set_SDF), "", "Signed-distance function to get the eff. dt at the interface")
			.add_method("set_vel_potential", static_cast<void (T::*)(SmartPtr<function_type>)>(&T::set_vel_potential), "", "Potential of the velocity")
			.add_method("prepare_for_SDF", static_cast<void (T::*)(SmartPtr<function_type>,SmartPtr<function_type>)>(&T::prepare_for_SDF), "old#new", "Prepare for computation of SDF")
			.add_method("set_verbose", &T::set_verbose)
			.add_method("set_dt", &T::set_dt)
			.add_method("get_dt", &T::get_dt)
			.add_method("set_time_control", &T::set_time_control)
			.add_method("set_time_control_off", &T::set_time_control_off)
			.add_method("set_time", &T::set_time)
			.add_method("get_time", &T::get_time)
			.add_method("get_CFL", &T::get_CFL)
			.add_method("set_divfree", &T::set_divfree)
			.add_method("set_first_order", &T::set_first_order)
			.add_method("set_antideriv_src", &T::set_antideriv_src)
			.add_method("set_elem_vel_vec", &T::set_elem_vel_vec)
			.add_method("set_nr_of_steps", &T::set_nr_of_steps)
			.add_method("set_delta", &T::set_delta)
			.add_method("set_gamma", &T::set_gamma)
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
			.add_method("set_source_neg", static_cast<void (T::*)(number)>(&T::set_source_neg), "", "Source (for LSF < 0)")
			.add_method("set_source_pos", static_cast<void (T::*)(number)>(&T::set_source_pos), "", "Source (for LSF > 0)")
			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity vector field")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity vector field")
#endif
			.add_method("set_normal_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_normal_velocity), "", "Normal velocity field")
#ifdef UG_FOR_LUA
			.add_method("set_normal_velocity", static_cast<void (T::*)(const char*)>(&T::set_normal_velocity), "", "Normal velocity field")
#endif
			.add_method("set_interface_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_interface_data), "", "Dirichlet BC at the interface")
#ifdef UG_FOR_LUA
			.add_method("set_interface_data", static_cast<void (T::*)(const char*)>(&T::set_interface_data), "", "Dirichlet BC at the interface")
#endif
			.add_method("compute_normal_vel", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >, SmartPtr<function_type>)>(&T::compute_normal_vel), "VelField#NormalVel", "Compute normal velocity")
#ifdef UG_FOR_LUA
			.add_method("compute_normal_vel", static_cast<void (T::*)(const char*, SmartPtr<function_type>)>(&T::compute_normal_vel), "VelField#NormalVel", "Compute normal velocity")
#endif
			.add_method("append_vertical_to_normal_vel", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, SmartPtr<function_type>)>(&T::append_vertical_to_normal_vel), "ZVelField#NormalVel", "Appends a vertical part to the normal velocity")
			.add_method("append_vertical_to_normal_vel", static_cast<void (T::*)(SmartPtr<function_type>, const char *, SmartPtr<function_type>)>(&T::append_vertical_to_normal_vel), "ZVelGF#cmp#NormalVel", "Appends a vertical part to the normal velocity")
			.add_method("set_outflow_boundary", &T::set_outflow_boundary)
			.add_method("set_dirichlet_boundary", &T::set_dirichlet_boundary)
			.add_method("set_dirichlet_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_dirichlet_data), "", "Dirichlet BC")
			.add_method("set_dirichlet_data", static_cast<void (T::*)(number)>(&T::set_dirichlet_data), "", "Dirichlet BC")
#ifdef UG_FOR_LUA
			.add_method("set_dirichlet_data", static_cast<void (T::*)(const char*)>(&T::set_dirichlet_data), "", "Dirichlet BC")
#endif
			.add_method("set_limiter", static_cast<void (T::*)(bool)>(&T::set_limiter))
			.add_method("compare_lsf_with", static_cast<void (T::*)(SmartPtr<function_type>, number)>(&T::compare_lsf_with))
			.add_method("save_CourantNumber_to", static_cast<void (T::*)(SmartPtr<function_type>)>(&T::save_CourantNumber_to))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HiResFluxBasedLSM", tag);
	}
	
//	LSVolume
	{
		typedef LSVolume<function_type> T;
		string name = string("LSVolume").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*) (SmartPtr<function_type>)> ("LSF")
			.add_method ("on_subsets", static_cast<void (T::*)(const char*)> (&T::on_subsets), "subsets", "Restrict to subsets")
			.add_method ("check_positivity", static_cast<void (T::*) (bool)> (&T::check_positivity), "flag", "Check positivity on/off")
			.add_method ("compute", static_cast<void (T::*)()> (&T::compute), "", "Compute the volumes")
			.add_method ("volume_plus", static_cast<number (T::*)() const> (&T::volume_plus), "", "Volume of the positive part")
			.add_method ("volume_plus_in_subsets", static_cast<number (T::*)(const char *) const> (&T::volume_plus_in_subsets), "subsets", "Volume of the positive part in subsets")
			.add_method ("volume_minus", static_cast<number (T::*)() const> (&T::volume_minus), "", "Volume of the negative part")
			.add_method ("volume_minus_in_subsets", static_cast<number (T::*)(const char *) const> (&T::volume_minus_in_subsets), "subsets", "Volume of the negative part in subsets")
			.add_method ("print_details", static_cast<void (T::*)() const> (&T::print_details), "", "Prints the detailed volumes")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LSVolume", tag);
	}
	
//	LSFbyRaster
	{
		typedef LSFbyRaster<function_type> T;
		string name = string("LSFbyRaster").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*) (const char*)> ("FileName")
			.add_method ("interpolate_to", static_cast<void (T::*)(SmartPtr<function_type>)> (&T::interpolate_to), "LSF", "Compute the LSF")
			.add_method ("set_relative_to", static_cast<void (T::*)(const char*)> (&T::set_relative_to), "subsets", "Consider raster values as relative to the top")
			.add_method ("set_relative_to", static_cast<void (T::*)(const char*,int)> (&T::set_relative_to), "subsets#gridLevel", "Consider raster values as relative to the top")
			.add_method ("set_rel_tolerance", static_cast<void (T::*)(number)> (&T::set_rel_tolerance), "tolerance", "Tolerance for the top ray tracer")
			.add_method ("set_rel_default", static_cast<void (T::*)(number)> (&T::set_rel_default), "value", "Default value for the z-coord. of the top")
			.add_method("set_local_top_faces_only", &T::set_local_top_faces_only, "enable", "In parallel computations, only local top faces are used for relative initalizations")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LSFbyRaster", tag);
	}

//	VelByNormalVel
	{
			string name = string("VelByNormalVel").append(suffix);
			typedef VelByNormalVel<function_type> T;
			typedef CplUserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<function_type>, SmartPtr<function_type>)>("NormalVel#LSF")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "VelByNormalVel", tag);
	}
	
//	EikonalVel
	{
			string name = string("EikonalVel").append(suffix);
			typedef EikonalVel<function_type> T;
			typedef CplUserData<MathVector<dim>, dim> TBase;
			
			reg.add_class_<T, TBase> (name, grp)
				.template add_constructor<void (*)(SmartPtr<function_type>)>("NormalVel#LSF")
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "EikonalVel", tag);
	}
	
} // end Domain Algebra

}; // end Functionality
} // end namespace LevelSet


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_LevelSet(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/LevelSet");
	typedef LevelSet::Functionality Functionality;

	try{
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
