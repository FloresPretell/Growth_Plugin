// created by Christian Wehner

// extern headers
#include <iostream>
#include <sstream>
#include <string>

#include "bridge/bridge.h"

#include "lib_algebra/cpu_algebra_types.h"

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#include <boost/function.hpp>
#include "level_set.h"

using namespace std;

namespace ug{
using namespace ug::bridge;

template <typename TDomain, typename TAlgebra>
static void Register__Algebra_Domain(bridge::Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, SurfaceDoFDistribution, TAlgebra> function_type;

//	group string
	stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgSuffix = bridge::GetDomainSuffix<TDomain>();
	dimAlgSuffix.append(GetAlgebraSuffix<TAlgebra>());

	string dimAlgTag = GetDomainTag<TDomain>();
	dimAlgTag.append(GetAlgebraTag<TAlgebra>());


// 	FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		typedef typename function_type::domain_type domain_type;
		typedef boost::function<void (number& value,
													  const MathVector<domain_type::dim>& x,
													  number time)> NumberFunctor;
		string name = string("FV1LevelSetDisc").append(dimAlgSuffix);
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
						.add_method("set_vel_x", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_vel_x))
						.add_method("set_vel_y", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_vel_y))
						.add_method("set_vel_z", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_vel_z))
						.add_method("set_source", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_source))
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
						.add_method("set_source", static_cast<void (T::*)()>(&T::set_source))
						.add_method("set_dirichlet_data", static_cast<void (T::*)(number)>(&T::set_dirichlet_data))
						.add_method("set_dirichlet_data", static_cast<void (T::*)(const NumberFunctor&)>(&T::set_dirichlet_data))
						.add_method("set_dirichlet_data", static_cast<void (T::*)()>(&T::set_dirichlet_data))
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
						.add_method("compute_error", &T::compute_error);
				reg.add_class_to_group(name, "FV1LevelSetDisc", dimAlgTag);
	}

}

template <typename TDomain>
static void RegisterIElemDiscs(bridge::Registry& reg, string grp)
{
//	dimension of domain
	static const int dim = TDomain::dim;

//	suffix and tag
	string dimSuffix = GetDomainSuffix<TDomain>();
	string dimTag = GetDomainTag<TDomain>();

// NOTHING TO REGISTER HERE
}


template <typename TAlgebra>
static bool Register__Algebra(bridge::Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
		Register__Algebra_Domain<Domain1d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_2
		Register__Algebra_Domain<Domain2d, TAlgebra>(reg, grp);
#endif
#ifdef UG_DIM_3
		Register__Algebra_Domain<Domain3d, TAlgebra>(reg, grp);
#endif

	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in Register__Algebra_DoFDistribution: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}


extern "C" void
InitUGPlugin_LevelSetPlugin(ug::bridge::Registry* reg, std::string parentGroup)
{
	std::string grp(parentGroup); grp.append("LevelSet/");

	bool bReturn = true;
#ifdef UG_CPU_1
	bReturn &= Register__Algebra<CPUAlgebra>(*reg, grp);
#endif
#ifdef UG_CPU_2
	bReturn &= Register__Algebra<CPUBlockAlgebra<2> >(*reg, grp);
#endif
#ifdef UG_CPU_3
	bReturn &= Register__Algebra<CPUBlockAlgebra<3> >(*reg, grp);
#endif
#ifdef UG_CPU_4
	bReturn &= Register__Algebra<CPUBlockAlgebra<4> >(*reg, grp);
#endif
#ifdef UG_CPU_VAR
	bReturn &= Register__Algebra<CPUVariableBlockAlgebra >(*reg, grp);
#endif

	try
	{
#ifdef UG_DIM_1
			RegisterIElemDiscs<Domain1d>(*reg, grp);
#endif
#ifdef UG_DIM_2
			RegisterIElemDiscs<Domain2d>(*reg, grp);
#endif
#ifdef UG_DIM_3
			RegisterIElemDiscs<Domain3d>(*reg, grp);
#endif
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibDisc_ElemDisc: "
				"Registration failed (using name " << ex.name << ").\n");
		return;
	}

	return;
}


}// end of namespace
