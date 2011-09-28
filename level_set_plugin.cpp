// created by Raphael Prohl

// extern headers
#include <iostream>
#include <sstream>
#include <string>

#include "bridge/ug_bridge.h"

#include "lib_algebra/algebra_selector.h"
#include "lib_algebra/algebra_types.h"

#include "lib_disc/dof_manager/conform/conform.h"
#include "lib_disc/dof_manager/p1conform/p1conform.h"

#include "lib_disc/spatial_discretization/constraints/constraint_interface.h"

#include <boost/function.hpp>
#include "level_set.h"

#ifdef UG_PARALLEL
	#include "lib_disc/parallelization/parallel_grid_function.h"
#endif

namespace ug{
using namespace ug::bridge;

template <typename TDomain, typename TAlgebra, typename TDoFDistribution>
void Register__Algebra_DoFDistribution_Domain(bridge::Registry& reg, string parentGroup)
{
//	typedef
	static const int dim = TDomain::dim;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;
	typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

#ifdef UG_PARALLEL
		typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
#else
		typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
#endif

//	group string
	stringstream grpSS; grpSS << parentGroup << "/" << dim << "d";
	string grp = grpSS.str();

//	suffix and tag
	string dimAlgDDSuffix = bridge::GetDomainSuffix<TDomain>();
	dimAlgDDSuffix.append(GetAlgebraSuffix<TAlgebra>());
	dimAlgDDSuffix.append(GetDoFDistributionSuffix<TDoFDistribution>());

	string dimAlgDDTag = GetDomainTag<TDomain>();
	dimAlgDDTag.append(GetAlgebraTag<TAlgebra>());
	dimAlgDDTag.append(GetDoFDistributionTag<TDoFDistribution>());


// 	FV1LevelSetDisc
	{
		typedef FV1LevelSetDisc<function_type> T;
		typedef typename function_type::domain_type domain_type;
		typedef boost::function<void (number& value,
													  const MathVector<domain_type::dim>& x,
													  number time)> NumberFunctor;
		string name = string("FV1LevelSetDisc").append(dimAlgDDSuffix);
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
				reg.add_class_to_group(name, "FV1LevelSetDisc", dimAlgDDTag);
	}

}

template <typename TDomain>
void RegisterIElemDiscs(bridge::Registry& reg, string grp)
{
//	dimension of domain
	static const int dim = TDomain::dim;

//	suffix and tag
	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

// NOTHING TO REGISTER HERE
}


template <typename TAlgebra, typename TDoFDistribution>
static bool Register__Algebra_DoFDistribution(bridge::Registry& reg, string parentGroup)
{
//	get group string
	string grp = parentGroup; grp.append("/Discretization");

	try
	{

#ifdef UG_DIM_1
//	Domain dependent part 1D
	{
		typedef Domain<1, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_2
//	Domain dependent part 2D
	{
		typedef Domain<2, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
#endif

#ifdef UG_DIM_3
//	Domain dependent part 3D
	{
		typedef Domain<3, MultiGrid, MGSubsetHandler> domain_type;
		Register__Algebra_DoFDistribution_Domain<domain_type, TAlgebra, TDoFDistribution>(reg, grp);
	}
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

template <typename TAlgebra>
static bool Register__Algebra(bridge::Registry& reg, string parentGroup)
{
	bool bReturn = true;
#ifdef DOF_P1
	bReturn &= Register__Algebra_DoFDistribution<TAlgebra, P1DoFDistribution>(reg, parentGroup);
#endif
#ifdef DOF_GEN
	bReturn &= Register__Algebra_DoFDistribution<TAlgebra, DoFDistribution >(reg, parentGroup);
#endif

	return bReturn;
}


extern "C" void InitUGPlugin(ug::bridge::Registry* reg)
{

	string parentGroup("/ug4/Plugins");

	bool bReturn = true;
	bReturn &= Register__Algebra<CPUAlgebra>(*reg, parentGroup);
//	bReturn &= Register__Algebra<CPUBlockAlgebra<2> >(*reg, parentGroup);
	bReturn &= Register__Algebra<CPUBlockAlgebra<3> >(*reg, parentGroup);
//	bReturn &= Register__Algebra<CPUBlockAlgebra<4> >(*reg, parentGroup);
//	bReturn &= Register__Algebra<CPUVariableBlockAlgebra >(*reg, parentGroup);

	try
	{
	//	get group string
		string grp("/ug4/Plugins/Discretization");

#ifdef UG_DIM_1
	//	Domain dependend part 1D
			RegisterIElemDiscs<Domain1d>(*reg, grp);
#endif

#ifdef UG_DIM_2
	//	Domain dependend part 2D
			RegisterIElemDiscs<Domain2d>(*reg, grp);
#endif

#ifdef UG_DIM_3
	//	Domain dependend part 3D
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
