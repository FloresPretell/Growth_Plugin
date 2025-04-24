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
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"
#include "lib_disc/spatial_disc/dom_disc_embb.h"

// plugin headers
#include "kappla_discr.h"
#include "Export.h"
//#include "Tensor_discr.h"

#include "ls_concentration_dependent_velocity_linker.h"
#include "ls_extend_velocity_linker.h"

#include "ls_tubulin_velocity_linker.h"
#include "ls_interface_influx_linker.h"
#include "ls_interface_influx_outside_linker.h"
#include "ls_interface_influx_top_based_curvature_linker.h"
#include "ls_interface_reaction_linker.h" // top based curvature also
#include "ls_interface_reaction_linker_inside_curv.h"							 
#include "ls_initial_value.h"
#include "ls_initial_value_interface.h"
#include "ls_tensor_linker.h"

using namespace std;
using namespace ug::bridge;

namespace ug
{
	namespace NeuroGrowth
	{

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
			static void Domain(Registry &reg, string grp)
			{
				static const int dim = TDomain::dim;
				typedef CPUAlgebra TLSFAlgebra;
				typedef GridFunction<TDomain, TLSFAlgebra> TLSFct;
				string suffix = GetDomainSuffix<TDomain>();
				string tag = GetDomainTag<TDomain>();

				// 	KAPPLA : calculate curvature
				{
					typedef Kappla_LS<TLSFct> T;
					typedef DebugWritingObject<TLSFAlgebra> TBase;

					// typedef typename function_type::domain_type domain_type;
					string name = string("Kappla_LS").append(suffix);
					reg.add_class_<T, TBase>(name, grp)
						 .add_constructor()

						 .add_method("set_solutions", static_cast<void (T::*)(SmartPtr<TLSFct>, SmartPtr<TLSFct>)>(&T::set_solutions), "old#new", "Solutions at the time levels")
						 .add_method("set_LSF", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_LSF), "", "Level-Set function to specity the interface")
						 .add_method("set_SDF", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_SDF), "", "Signed-distance function to get the eff. dt at the interface")
						 .add_method("set_vel_potential", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_vel_potential), "", "Potential of the velocity")

						 .add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim>>)>(&T::set_velocity), "", "Velocity vector field")
						 .add_method("set_normal_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_normal_velocity), "", "Normal velocity field")
						 .add_method("set_interface_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_interface_data), "", "Dirichlet BC at the interface")
						 .add_method("compute_kappla", static_cast<void (T::*)(SmartPtr<TLSFct>, SmartPtr<TLSFct>)>(&T::compute_kappla), "Gridfunction of level set", "Gridfunction to save Kappla")

						 //.add_method("set_outflow_boundary", &T::set_outflow_boundary)

						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "Kappla_LS", tag);
				}

				/*
				// 	KAPPLA : calculate curvature
				{
					typedef TensorDiffusion_LS<TLSFct> T;
					typedef DebugWritingObject<TLSFAlgebra> TBase;

					// typedef typename function_type::domain_type domain_type;
					string name = string("TensorDiffusion_LS").append(suffix);
					reg.add_class_<T, TBase>(name, grp)
						 .add_constructor()

						 .add_method("set_solutions", static_cast<void (T::*)(SmartPtr<TLSFct>, SmartPtr<TLSFct>)>(&T::set_solutions), "old#new", "Solutions at the time levels")
						 .add_method("set_LSF", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_LSF), "", "Level-Set function to specity the interface")
						 .add_method("set_SDF", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_SDF), "", "Signed-distance function to get the eff. dt at the interface")
						 .add_method("set_vel_potential", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_vel_potential), "", "Potential of the velocity")

						 .add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim>>)>(&T::set_velocity), "", "Velocity vector field")
						 .add_method("set_normal_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_normal_velocity), "", "Normal velocity field")
						 .add_method("set_interface_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_interface_data), "", "Dirichlet BC at the interface")
						 .add_method("ComputeDiffusionTensorBeltrami", static_cast<void (T::*)(SmartPtr<TLSFct>, SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >)>(&T::ComputeDiffusionTensorBeltrami), "Gridfunction of level set", "Gridfunction to save Kappla")

						 //.add_method("set_outflow_boundary", &T::set_outflow_boundary)

						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "TensorDiffusion_LS", tag);
				}
				*/
				// 	Export : Export a linker as a gridfunction
				{
					typedef Export<TLSFct> T;
					typedef DebugWritingObject<TLSFAlgebra> TBase;

					// typedef typename function_type::domain_type domain_type;
					string name = string("Export").append(suffix);
					reg.add_class_<T, TBase>(name, grp)
						 .add_constructor()

						 .add_method("set_solutions", static_cast<void (T::*)(SmartPtr<TLSFct>, SmartPtr<TLSFct>)>(&T::set_solutions), "old#new", "Solutions at the time levels")
						 .add_method("set_LSF", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_LSF), "", "Level-Set function to specity the interface")
						 .add_method("set_SDF", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_SDF), "", "Signed-distance function to get the eff. dt at the interface")
						 .add_method("set_vel_potential", static_cast<void (T::*)(SmartPtr<TLSFct>)>(&T::set_vel_potential), "", "Potential of the velocity")
						 .add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim>>)>(&T::set_velocity), "", "Velocity vector field")
						 .add_method("set_normal_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_normal_velocity), "", "Normal velocity field")
						 .add_method("set_interface_data", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_interface_data), "", "Dirichlet BC at the interface")
						 .add_method("Save_Linker_to_Gridfuntion", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim>>, SmartPtr<TLSFct>)>(&T::Save_Linker_to_Gridfuntion), "VelField#NormalVel", "Compute normal velocity")
						 .add_method("Save_Linker_to_Gridfuntion_by_component", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim>>, SmartPtr<TLSFct>, const char* )>(&T::Save_Linker_to_Gridfuntion_by_component), "VelField#NormalVel", "Compute normal velocity","Componeten a elegir")
						 
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "Export", tag);
				}


				
			}




			/**
			* Function called for the registration of Dimension dependent parts.
			* All Functions and Classes depending on the Dimension
			* are to be placed here when registering. The method is called for all
			* available Dimension types, based on the current build options.
			*
			* @param reg				registry
			* @param parentGroup		group for sorting of functionality
			*/
			template <int dim>
			static void Dimension(Registry& reg, string grp)
			{
				string suffix = GetDimensionSuffix<dim>();
				string tag = GetDimensionTag<dim>();

				//	Linker for TENSOR
				{
					string name = string("LSTensorLinker").append(suffix);
					typedef LSTensorLinker<dim> T;
					typedef DependentUserData<MathMatrix<dim, dim>, dim> TBase;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_gradient_stationary_difussion", &T::set_gradient_stationary_difussion)
						 .add_method("set_diffusion_coeff", static_cast<void (T::*)(number)>(&T::set_diffusion_coeff))
						 .add_constructor()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSTensorLinker", tag);
				}


			} // end Domain Algebra


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
			static void DomainAlgebra(Registry &reg, string grp)
			{
				string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
				string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

				// typedef typename TAlgebra::vector_type vector_type;
				// typedef typename TAlgebra::matrix_type matrix_type;
				// typedef ApproximationSpace<TDomain> approximation_space_type;
				typedef GridFunction<TDomain, TAlgebra> function_type;
				static const int dim = TDomain::dim;

				//	Special linker for the Darcy velocity
				{
					string name = string("LSConcentrationDepentVelocity").append(suffix);
					typedef LSConcentrationDepentVelocity<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_calcium", static_cast<void (T::*)(number)>(&T::set_calcium))
						 .add_method("set_calcium", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_calcium))
						 .add_method("set_tubuline", static_cast<void (T::*)(number)>(&T::set_tubuline))
						 .add_method("set_tubuline", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_tubuline))
						 .add_method("set_MAPu", static_cast<void (T::*)(number)>(&T::set_MAPu))
						 .add_method("set_MAPu", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_MAPu))
						 .add_method("set_MAPb", static_cast<void (T::*)(number)>(&T::set_MAPb))
						 .add_method("set_MAPb", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_MAPb))
						 .add_method("set_MAPp", static_cast<void (T::*)(number)>(&T::set_MAPp))
						 .add_method("set_MAPp", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_MAPp))
						 .add_method("set_Inhibition", static_cast<void (T::*)(number)>(&T::set_Inhibition))
						 .add_method("set_Inhibition", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Inhibition))
						 .add_method("set_LevelSet_gradient", &T::set_LevelSet_gradient)
						 .add_method("set_gradient_stationary_difussion", &T::set_gradient_stationary_difussion)
						 .add_method("set_Inhibition_gradient", &T::set_Inhibition_gradient)
						 .add_method("set_Curvature", static_cast<void (T::*)(number)>(&T::set_Curvature))
						 .add_method("set_Curvature", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Curvature))
						 .add_method("set_interval_min", static_cast<void (T::*)(number)>(&T::set_interval_min))
						 .add_method("set_interval_max", static_cast<void (T::*)(number)>(&T::set_interval_max))
						 .add_method("set_interval_min_Calcium", static_cast<void (T::*)(number)>(&T::set_interval_min_Calcium))
						 .add_method("set_interval_min_Inhibition", static_cast<void (T::*)(number)>(&T::set_interval_min_Inhibition))
						 .add_method("set_interval_min_Inhibition_sign", static_cast<void (T::*)(number)>(&T::set_interval_min_Inhibition_sign))
						 .add_method("set_magnitud_velocity", static_cast<void (T::*)(number)>(&T::set_magnitud_velocity))
						 .template add_constructor<void (*)(SmartPtr<TExtrapol>)>("DomainDisc")
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSConcentrationDepentVelocity", tag);
				}

				//	Special linker for calculate the velocity of the extencion of the calculated velocity of the interface
				{
					string name = string("LSExtendVelocity").append(suffix);
					typedef LSExtendVelocity<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_gradient", &T::set_gradient)
						 .add_method("set_direction", &T::set_direction)
						 .template add_constructor<void (*)(SmartPtr<TExtrapol>)>("DomainDisc")
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSExtendVelocity", tag);
				}


				//	Special linker for initial values
				{
					string name = string("LSInitialValue").append(suffix);
					typedef LSInitialValue<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_initial_value", static_cast<void (T::*)(number)>(&T::set_initial_value))
						 .template add_constructor<void (*)(SmartPtr<TExtrapol>)>("DomainDisc")
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSInitialValue", tag);
				}

				
				//	Special linker for initial values IN THE INTERFAXE
				{
					string name = string("LSInitialValueInterface").append(suffix);
					typedef LSInitialValueInterface<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_initial_value", static_cast<void (T::*)(number)>(&T::set_initial_value))
						 .template add_constructor<void (*)(SmartPtr<TExtrapol>)>("DomainDisc")
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSInitialValueInterface", tag);
				}


				//	Linker for obtain the direction of the tubulin
				{
					string name = string("LSTubulinVelocity").append(suffix);
					typedef LSTubulinVelocity<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_gradient_stationary_difussion", &T::set_gradient_stationary_difussion)
						 .add_method("set_domain_discretizacion", &T::set_domain_discretizacion)
						 .add_method("set_magnitud", static_cast<void (T::*)(number)>(&T::set_magnitud))
						 .template add_constructor<void (*)()>()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSTubulinVelocity", tag);
				}



				//	Linker for give the value to influx for the boundary (cut elements)
				{
					string name = string("LSInterfaceInflux").append(suffix);
					typedef LSInterfaceInflux<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_LevelSet_gradient", &T::set_LevelSet_gradient)
						 .add_method("set_Tubuline", static_cast<void (T::*)(number)>(&T::set_Tubuline))
						 .add_method("set_Tubuline", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Tubuline))
						 .add_method("set_domain_discretizacion", &T::set_domain_discretizacion)
						 .add_method("set_magnitud_influx", static_cast<void (T::*)(number)>(&T::set_magnitud_influx))
						 //.add_method("set_domain_disc_1d"       , &T::set_domain_disc_1d, "", "domainDisc","Set the 1d cable domain discretization.")
						 .template add_constructor<void (*)()>()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSInterfaceInflux", tag);
				}

				//	Linker for give the value to influx for the boundary (cut elements) to outside the domain
				{
					string name = string("LSInterfaceInfluxoutside").append(suffix);
					typedef LSInterfaceInfluxoutside<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_LevelSet_gradient", &T::set_LevelSet_gradient)
						 .add_method("set_Tubuline", static_cast<void (T::*)(number)>(&T::set_Tubuline))
						 .add_method("set_Tubuline", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Tubuline))
 						 .add_method("set_Curvature", static_cast<void (T::*)(number)>(&T::set_Curvature))
						 .add_method("set_Curvature", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Curvature))
						 .add_method("set_domain_discretizacion", &T::set_domain_discretizacion)
						 .add_method("set_carrying_capacity", static_cast<void (T::*)(number)>(&T::set_carrying_capacity))
						 .add_method("set_growth_rate", static_cast<void (T::*)(number)>(&T::set_growth_rate))
						 .add_method("set_interval_min", static_cast<void (T::*)(number)>(&T::set_interval_min))
						 .add_method("set_interval_max", static_cast<void (T::*)(number)>(&T::set_interval_max))
							//.add_method("set_domain_disc_1d"       , &T::set_domain_disc_1d, "", "domainDisc","Set the 1d cable domain discretization.")
						 .template add_constructor<void (*)()>()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSInterfaceInfluxoutside", tag);
				}

				//	Linker for give the value to influx for the boundary (cut elements) - This will select the top of the branch with a curvature specify
				{
					string name = string("LSInfluxTopBranch").append(suffix);
					typedef LSInfluxTopBranch<TDomain, TAlgebra> T;
					typedef DependentUserData<MathVector<dim>, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_LevelSet_gradient", &T::set_LevelSet_gradient)
						 .add_method("set_magnitud_influx", static_cast<void (T::*)(number)>(&T::set_magnitud_influx))	 
						 .add_method("set_Curvature", static_cast<void (T::*)(number)>(&T::set_Curvature))
						 .add_method("set_Curvature", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Curvature))
						 .add_method("set_domain_discretizacion", &T::set_domain_discretizacion)
						 .add_method("set_interval_min", static_cast<void (T::*)(number)>(&T::set_interval_min))
						 .add_method("set_interval_max", static_cast<void (T::*)(number)>(&T::set_interval_max))
						 //.add_method("set_domain_disc_1d"       , &T::set_domain_disc_1d, "", "domainDisc","Set the 1d cable domain discretization.")
						 .template add_constructor<void (*)()>()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSInfluxTopBranch", tag);
				}



				//	Linker for give the value to REACTION for the boundary (cut elements) - This will select the top of the branch with a curvature specify
				{
					string name = string("LSReactionTopBranch").append(suffix);
					typedef LSReactionTopBranch<TDomain, TAlgebra> T;
					typedef DependentUserData<number, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_FirstTerm", static_cast<void (T::*)(number)>(&T::set_FirstTerm))
						 .add_method("set_FirstTerm", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_FirstTerm))
						 .add_method("set_SecondTerm", static_cast<void (T::*)(number)>(&T::set_SecondTerm))
						 .add_method("set_SecondTerm", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_SecondTerm))
						 .add_method("set_LevelSet_gradient", &T::set_LevelSet_gradient)
						 .add_method("set_Curvature", static_cast<void (T::*)(number)>(&T::set_Curvature))
						 .add_method("set_Curvature", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Curvature))
						 .add_method("set_velocity_deformation", static_cast<void (T::*)(number)>(&T::set_velocity_deformation))
						 .add_method("set_velocity_deformation", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_velocity_deformation))
						 .add_method("set_domain_discretizacion", &T::set_domain_discretizacion)
						 .add_method("set_interval_min", static_cast<void (T::*)(number)>(&T::set_interval_min))
						 .add_method("set_interval_max", static_cast<void (T::*)(number)>(&T::set_interval_max))
						 .add_method("set_magnitud_influx", static_cast<void (T::*)(number)>(&T::set_magnitud_influx))
						 .add_method("set_constant_div", static_cast<void (T::*)(number)>(&T::set_constant_div))
						 .add_method("set_constant_addition", static_cast<void (T::*)(number)>(&T::set_constant_addition))

						 //.add_method("set_domain_disc_1d"       , &T::set_domain_disc_1d, "", "domainDisc","Set the 1d cable domain discretization.")
						 .template add_constructor<void (*)()>()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSReactionTopBranch", tag);
				}

				//	Linker for give the value to REACTION for the boundary (cut elements) - This will select the top of the branch with a curvature specify
				{
					string name = string("LSReactionInsideCurv").append(suffix);
					typedef LSReactionInsideCurv<TDomain, TAlgebra> T;
					typedef DependentUserData<number, dim> TBase;
					typedef IInterfaceExtrapolation<TDomain, TAlgebra> TExtrapol;
					reg.add_class_<T, TBase>(name, grp)
						 .add_method("set_FirstTerm", static_cast<void (T::*)(number)>(&T::set_FirstTerm))
						 .add_method("set_FirstTerm", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_FirstTerm))
						 .add_method("set_SecondTerm", static_cast<void (T::*)(number)>(&T::set_SecondTerm))
						 .add_method("set_SecondTerm", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_SecondTerm))
						 .add_method("set_Curvature", static_cast<void (T::*)(number)>(&T::set_Curvature))
						 .add_method("set_Curvature", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_Curvature))
						 .add_method("set_velocity_deformation", static_cast<void (T::*)(number)>(&T::set_velocity_deformation))
						 .add_method("set_velocity_deformation", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim>>)>(&T::set_velocity_deformation))
						 .add_method("set_domain_discretizacion", &T::set_domain_discretizacion)
						 .add_method("set_interval_min", static_cast<void (T::*)(number)>(&T::set_interval_min))
						 .add_method("set_interval_max", static_cast<void (T::*)(number)>(&T::set_interval_max))
						 .add_method("set_magnitud_influx", static_cast<void (T::*)(number)>(&T::set_magnitud_influx))
						 .add_method("set_constant_div", static_cast<void (T::*)(number)>(&T::set_constant_div))
						 .add_method("set_constant_addition", static_cast<void (T::*)(number)>(&T::set_constant_addition))

						 //.add_method("set_domain_disc_1d"       , &T::set_domain_disc_1d, "", "domainDisc","Set the 1d cable domain discretization.")
						 .template add_constructor<void (*)()>()
						 .set_construct_as_smart_pointer(true);
					reg.add_class_to_group(name, "LSReactionInsideCurv", tag);
				}
			}



		}; // end Functionality
	} // end namespace NeuroGrowth

	/**
	 * This function is called when the plugin is loaded.
	 */
	extern "C" void
	InitUGPlugin_NeuroGrowth(Registry *reg, string grp)
	{
		grp.append("/SpatialDisc/NeuroGrowth");
		typedef NeuroGrowth::Functionality Functionality;

		try
		{
			RegisterDimensionDependent<Functionality>(*reg,grp);
			RegisterDomain2d3dDependent<Functionality>(*reg, grp);
			RegisterDomain2d3dAlgebraDependent<Functionality>(*reg, grp);
		}
		UG_REGISTRY_CATCH_THROW(grp);
	}

} // namespace ug
