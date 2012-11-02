// created by Christian Wehner
/**
 * LevelSet plugin
 *
 * created by Christian Wehner
 */

#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "level_set.h"

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
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> function_type;
	static const int dim = TDomain::dim;

// 	FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		typedef typename function_type::domain_type domain_type;
		string name = string("FV1LevelSetDisc").append(suffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("set_dt", &T::set_dt)
			.add_method("set_vel_scale", &T::set_vel_scale)
			.add_method("set_reinit", &T::set_reinit)
			.add_method("set_divfree_bool", &T::set_divfree_bool)
			.add_method("set_info", &T::set_info)
			.add_method("set_timestep_nr",&T::set_timestep_nr)
			.add_method("set_nr_of_steps", &T::set_nr_of_steps)
			.add_method("advect_lsf", &T::advect_lsf)
			.add_method("add_post_process", &T::add_post_process)
			.add_method("set_delta", &T::set_delta)
			.add_method("set_gamma", &T::set_gamma)
			.add_method("set_limiter",&T::set_limiter)
			.add_method("set_dirichlet_boundary",&T::set_dirichlet_boundary)
			.add_method("set_outflow_boundary",&T::set_outflow_boundary)
			.add_method("init_function", &T::init_function)
			.add_method("set_source", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif
			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(number)>(&T::set_velocity), "", "Vel_x")
			.add_method("set_velocity", static_cast<void (T::*)(number,number)>(&T::set_velocity), "", "Vel_x, Vel_y")
			.add_method("set_velocity", static_cast<void (T::*)(number,number,number)>(&T::set_velocity), "", "Vel_x, Vel_y, Vel_z")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
#endif
/*			.add_method("set_vel_x", static_cast<void (T::*)(SmartPtr<UserData<number,dim> >)>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)(SmartPtr<UserData<number,dim> >)>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)(SmartPtr<UserData<number,dim> >)>(&T::set_vel_z))
			.add_method("set_source", static_cast<void (T::*)(SmartPtr<UserData<number,dim> >)>(&T::set_source))
			.add_method("set_vel_x", static_cast<void (T::*)(function_type&)>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)(function_type&)>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)(function_type&)>(&T::set_vel_z))
			.add_method("set_vel_x", static_cast<void (T::*)(number)>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)(number)>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)(number)>(&T::set_vel_z))
			.add_method("set_vel_x", static_cast<void (T::*)()>(&T::set_vel_x))
			.add_method("set_vel_y", static_cast<void (T::*)()>(&T::set_vel_y))
			.add_method("set_vel_z", static_cast<void (T::*)()>(&T::set_vel_z))
			.add_method("set_source", static_cast<void (T::*)(function_type&)>(&T::set_source))
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source))
			.add_method("set_source", static_cast<void (T::*)()>(&T::set_source))*/
			.add_method("set_dirichlet_data", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_dirichlet_data), "", "Source")
			.add_method("set_dirichlet_data", static_cast<void (T::*)(number)>(&T::set_dirichlet_data), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_dirichlet_data", static_cast<void (T::*)(const char*)>(&T::set_dirichlet_data), "", "Source")
#endif
			.add_method("compute_normal",&T::compute_normal)
			.add_method("compute_dnormal",&T::compute_dnormal)
			.add_method("compute_ddnormal",&T::compute_ddnormal)
			.add_method("fill_v_vec",&T::fill_v_vec)
			.add_method("get_time",&T::get_time)
			.add_method("runtimetest", &T::runtimetest)
			.add_method("init_ls_subsets",&T::init_ls_subsets)
			.add_method("create_ls_subsets",&T::create_ls_subsets)
			.add_method("update_ls_subsets",&T::update_ls_subsets)
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
			.add_method("compute_curvature", &T::computeElementCurvatureOnGrid2d);
		reg.add_class_to_group(name, "FV1LevelSetDisc", tag);
	}

}

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
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
