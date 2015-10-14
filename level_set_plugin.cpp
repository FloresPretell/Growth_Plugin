/**
 * LevelSet plugin
 *
 * created by Christian Wehner
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

	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> function_type;
	static const int dim = TDomain::dim;

// 	FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		typedef typename function_type::domain_type domain_type;
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
		typedef typename function_type::domain_type domain_type;
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
		typedef typename function_type::domain_type domain_type;
		string name = string("HiResFluxBasedLSM").append(suffix);
		reg.add_class_<T>(name, grp)
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
			.add_method("append_to_normal_vel", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >, SmartPtr<function_type>)>(&T::append_to_normal_vel), "NormVelField#NormalVel", "Appends a horizontal part to the normal velocity")
			.add_method("append_to_normal_vel", static_cast<void (T::*)(SmartPtr<function_type>, const char *, SmartPtr<function_type>)>(&T::append_to_normal_vel), "NormVelGF#cmp#NormalVel", "Appends a horizontal part to the normal velocity")
			.add_method("set_outflow_boundary", &T::set_outflow_boundary)
			.add_method("set_dirichlet_boundary", &T::set_dirichlet_boundary)
			.add_method("set_dirichlet_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_dirichlet_data), "", "Dirichlet BC")
			.add_method("set_dirichlet_data", static_cast<void (T::*)(number)>(&T::set_dirichlet_data), "", "Dirichlet BC")
#ifdef UG_FOR_LUA
			.add_method("set_dirichlet_data", static_cast<void (T::*)(const char*)>(&T::set_dirichlet_data), "", "Dirichlet BC")
#endif
			.add_method("set_limiter", static_cast<void (T::*)(bool)>(&T::set_limiter))
			.add_method("compare_lsf_with", static_cast<void (T::*)(SmartPtr<function_type>, number)>(&T::compare_lsf_with))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "HiResFluxBasedLSM", tag);
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
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
