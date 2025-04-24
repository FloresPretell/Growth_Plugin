/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
 * Based on the module by Andreas Vogel.
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

#ifndef __H__UG__PLUGINS__LEVEL_SET_EXTEND_VELOCITY_LINKER__
#define __H__UG__PLUGINS__LEVEL_SET_EXTEND_VELOCITY_LINKER__

// ug4 headers
#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/dom_disc_embb.h"


#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#endif


namespace ug
{
	namespace NeuroGrowth
	{



		template <typename TDomain, typename TAlgebra>
		class LSExtendVelocity
			 : public StdDataLinker<LSExtendVelocity<TDomain, TAlgebra>, MathVector<TDomain::dim>, TDomain::dim>
		{
			//	domain type
			typedef TDomain domain_type;

			//	algebra type
			typedef TAlgebra algebra_type;

			//	world dimension
			static const int dim = domain_type::dim;

			//	Base class type
			typedef StdDataLinker<LSExtendVelocity<domain_type, algebra_type>, MathVector<dim>, dim> base_type;

			//	extrapolation type
			typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

		public:
			LSExtendVelocity(SmartPtr<extrapol_type> spExtrapol) : 
				m_spExtrapolation(spExtrapol),
                m_spLevelSetGrad(NULL) , m_spDLevelSetGrad(NULL),
				m_spInhibitionGrad(NULL), m_spDInhibitionGrad(NULL) // poner la gradiente de la inhibicion
			{
				//	this linker needs exactly five input
				this->set_num_input(2);
			}

		private:
			//	checks whether the element is cut
			int elem_cut(
				 GridObject *elem,							 ///< the element to process
				 int si,											 ///< subset of the element
				 const MathVector<dim> vCornerCoords[], ///< coordinates of the corners
				 number time									 ///< the phisical time
			) const
			{
				if (m_spExtrapolation.valid())
				{
					const ReferenceObjectID roid = elem->reference_object_id();
					const DimReferenceElement<dim> &rRefElem = ReferenceElementProvider::get<dim>(roid);
					return ((extrapol_type *)m_spExtrapolation.get())->check_elem_lsf(rRefElem.num(0), elem, si, false, vCornerCoords, time);
				}
				return 1;
			}

		public:
			inline void evaluate(MathVector<dim> &value,
										const MathVector<dim> &globIP,
										number time, int si) const
			{
				UG_THROW("LSExtendVelocity: Element is necessary for the evaluation.");
			}

			template <int refDim>
			inline void evaluate(MathVector<dim> vValue[],
										const MathVector<dim> vGlobIP[],
										number time, int si,
										GridObject *elem,
										const MathVector<dim> vCornerCoords[],
										const MathVector<refDim> vLocIP[],
										const size_t nip,
										LocalVector *u,
										const MathMatrix<refDim, dim> *vJT = NULL) const
			{
				if (elem_cut(elem, si, vCornerCoords, time) > 0) // if under the free surface
				{
					for (size_t ip = 0; ip < nip; ++ip)
						vValue[ip] = 0.0;
					return;
				}


           		std::vector<MathVector<dim>> vInhibitionGrad(nip);
                std::vector<MathVector<dim>> vLevelSetGrad(nip);

            	(*m_spInhibitionGrad)(&vInhibitionGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
                (*m_spLevelSetGrad)(&vLevelSetGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);

				for (size_t ip = 0; ip < nip; ++ip)
				{
                    //	Calculate Magnitud
                    number magnitud = VecLength(vInhibitionGrad[ip]);                     // obtain the norm of the gradient of the LevelSet

                    //	Calculate Direccion = normal
                    MathVector<dim> normal;                                              // vector to save direction of the velocity
                    number gradientnorm = VecLength(vLevelSetGrad[ip]);                     // obtain the norm of the gradient of the LevelSet
                    VecScale(normal, vLevelSetGrad[ip], 1.0 / (gradientnorm + 0.00001)); // normalize the gradient of the LevelSet

                    MathVector<dim> VelExt;
                    VecScale(VelExt, normal, magnitud); // normalize the gradient of the LevelSet

                    VecScale(vValue[ip], VelExt, 1.0);
				}
			}

			template <int refDim>
			void eval_and_deriv(MathVector<dim> vDarcyVel[],
									  const MathVector<dim> vGlobIP[],
									  number time, int si,
									  GridObject *elem,
									  const MathVector<dim> vCornerCoords[],
									  const MathVector<refDim> vLocIP[],
									  const size_t nip,
									  LocalVector *u,
									  bool bDeriv,
									  int s,
									  std::vector<std::vector<MathVector<dim>>> vvvDeriv[],
									  const MathMatrix<refDim, dim> *vJT = NULL) const
			{
				if (elem_cut(elem, si, vCornerCoords, time) > 0) // if under the free surface
				{
					for (size_t ip = 0; ip < nip; ++ip)
						vDarcyVel[ip] = 0.0;
					if (!bDeriv || this->zero_derivative())
						return;
					this->set_zero(vvvDeriv, nip);
					return;
				}

				//	get the data of the ip series
				const MathVector<dim>* vInhibitionGrad = m_spInhibitionGrad->values(s);
                const MathVector<dim>* vLevelSetGrad = m_spLevelSetGrad->values(s);
				
                for (size_t ip = 0; ip < nip; ++ip)
				{
                    //	Calculate Magnitud
                    number magnitud = VecLength(vInhibitionGrad[ip]);                     // obtain the norm of the gradient of the LevelSet

                    //	Calculate Direccion = normal
                    MathVector<dim> normal;                                              // vector to save direction of the velocity
                    number gradientnorm = VecLength(vLevelSetGrad[ip]);                     // obtain the norm of the gradient of the LevelSet
                    VecScale(normal, vLevelSetGrad[ip], 1.0 / (gradientnorm + 0.00001)); // normalize the gradient of the LevelSet

                    MathVector<dim> VelExt;
                    VecScale(VelExt, normal, magnitud); // normalize the gradient of the LevelSet

                    VecScale(vDarcyVel[ip], VelExt, 1.0);
				}

				//	check if something to do
					if(!bDeriv || this->zero_derivative()) return;

				//	clear all derivative values
					this->set_zero(vvvDeriv, nip);
			}

		public:

			///	set LevelSet gradient import
			void set_gradient(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
			{
				m_spInhibitionGrad = data;
				m_spDInhibitionGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim>>();
				base_type::set_input(_LG_, data, data);
			}

			///	set the direction : this will use the gradient of a imagina
			void set_direction(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
			{
				m_spLevelSetGrad = data;
				m_spDLevelSetGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim>>();
				base_type::set_input(_DIR_, data, data);
			}

		protected:
			///	extrapolation by the level-set function
			SmartPtr<extrapol_type> m_spExtrapolation;

			///	import for LevelSet gradient
			static const size_t _LG_ = 0;
			SmartPtr<CplUserData<MathVector<dim>, dim>> m_spInhibitionGrad;
			SmartPtr<DependentUserData<MathVector<dim>, dim>> m_spDInhibitionGrad;

			///	import the gradient of a imaginary molecule = direccion
			static const size_t _DIR_ = 1;
			SmartPtr<CplUserData<MathVector<dim>, dim>> m_spLevelSetGrad;
			SmartPtr<DependentUserData<MathVector<dim>, dim>> m_spDLevelSetGrad;

		};

	} // end namespace NeuroGrowth
} // end namespace ug

#endif /* __H__UG__PLUGINS__LEVEL_SET_DARCY_VELOCITY_LINKER__ */
