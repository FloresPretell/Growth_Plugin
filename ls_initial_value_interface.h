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

#ifndef __H__UG__PLUGINS__LEVEL_SET_CONCENTRATION_INTERFACE_LINKER2__
#define __H__UG__PLUGINS__LEVEL_SET_CONCENTRATION_INTERFACE_LINKER2__

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
		class LSInitialValueInterface
			 : public StdDataLinker<LSInitialValueInterface<TDomain, TAlgebra>, MathVector<TDomain::dim>, TDomain::dim>
		{
			//	domain type
			typedef TDomain domain_type;

			//	algebra type
			typedef TAlgebra algebra_type;

			//	world dimension
			static const int dim = domain_type::dim;

			//	Base class type
			typedef StdDataLinker<LSInitialValueInterface<domain_type, algebra_type>, MathVector<dim>, dim> base_type;

			//	extrapolation type
			typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

		public:
			LSInitialValueInterface(SmartPtr<extrapol_type> spExtrapol) : 
				m_spExtrapolation(spExtrapol),
				m_spvalue(NULL)

			{
				//	this linker needs exactly five input
				this->set_num_input(1);
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
				UG_THROW("LSInitialValueInterface: Element is necessary for the evaluation.");
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
				if (elem_cut(elem, si, vCornerCoords, time) < 0 || elem_cut(elem, si, vCornerCoords, time) > 0)
				{
					for (size_t ip = 0; ip < nip; ++ip)
						vValue[ip] = 0.0;
					return;
				}

				// el problema es el mapeo dado que tienen diferente tamaño
				// vCalcium pide 3 -- se verifico que vGlobIP
				// m_spCalcium solo te entrega 2

				const number initialvalue = m_spvalue; //0.5

				for (size_t ip = 0; ip < nip; ++ip)
				{
						MathVector<dim> Vel;
						VecSet(Vel, initialvalue);
						VecScale(vValue[ip], Vel, 1.0);
					
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
				if (elem_cut(elem, si, vCornerCoords, time) < 0 || elem_cut(elem, si, vCornerCoords, time) > 0)
				{
					for (size_t ip = 0; ip < nip; ++ip)
						vDarcyVel[ip] = 0.0;
					if (!bDeriv || this->zero_derivative())
						return;
					this->set_zero(vvvDeriv, nip);
					return;
				}

				//	get the data of the ip series
				const float initialvalue = m_spvalue; //0.5

				for (size_t ip = 0; ip < nip; ++ip)
				{
						MathVector<dim> Vel;
						VecSet(Vel, initialvalue);
						VecScale(vDarcyVel[ip], Vel, 1.0);
				}

				//	check if something to do
					if(!bDeriv || this->zero_derivative()) return;

				//	clear all derivative values
					this->set_zero(vvvDeriv, nip);
			}

		public:
			///	set initial value : this vill be to generate the linker
			void set_initial_value(number data)
			{
				m_spvalue = data;
			}



		protected:
			///	extrapolation by the level-set function
			SmartPtr<extrapol_type> m_spExtrapolation;

			///	set initial value : this vill be to generate the linker		
            number m_spvalue;
		};

	} // end namespace NeuroGrowth
} // end namespace ug

#endif /* __H__UG__PLUGINS__LEVEL_SET_DARCY_VELOCITY_LINKER__ */
