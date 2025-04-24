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

#ifndef __H__UG__PLUGINS__LEVEL_SET_TUBULINE_VELOCITY_LINKER__
#define __H__UG__PLUGINS__LEVEL_SET_TUBULINE_VELOCITY_LINKER__

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
      class LSTubulinVelocity
          : public StdDataLinker<LSTubulinVelocity<TDomain, TAlgebra>, MathVector<TDomain::dim>, TDomain::dim>
      {
         //	domain type
         typedef TDomain domain_type;

         //	algebra type
         typedef TAlgebra algebra_type;

         //	world dimension
         static const int dim = domain_type::dim;

         //	Base class type
         typedef StdDataLinker<LSTubulinVelocity<domain_type, algebra_type>, MathVector<dim>, dim> base_type;

         //	extrapolation type
         typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

      public:
         LSTubulinVelocity() : m_spExtrapolation(NULL),
                               m_spMagnitud(NULL),
                               m_spDirection(NULL) 
         {
            //	this linker needs exactly five input
            this->set_num_input(0);
         }

      private:
         //	checks whether the element is cut
         int elem_cut(
             GridObject *elem,                      ///< the element to process
             int si,                                ///< subset of the element
             const MathVector<dim> vCornerCoords[], ///< coordinates of the corners
             number time                            ///< the phisical time
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
            UG_THROW("LSDarcyVelocityLinker: Element is necessary for the evaluation.");
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

            if (elem_cut(elem, si, vCornerCoords, time) < 0 ) //  outside (<0) 
            {
               for (size_t ip = 0; ip < nip; ++ip)
               {
                  MathVector<dim> Vel;
                  VecSet(Vel, 0);
                  VecScale(vValue[ip], Vel, 1.0);
               }
               return; // salir de de evaluate
            }

            std::vector<MathVector<dim>> vDirection(nip);

            (*m_spDirection)(&vDirection[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT);

				const number Magnitud = m_spMagnitud; // 0.0005 // cantidad de concentración que entrará por las paredes

            for (size_t ip = 0; ip < nip; ++ip)
            {
               //	Calculate Direccion
               MathVector<dim> direccion;                                              // vector to save direction of the velocity
               number gradientnorm = VecLength(vDirection[ip]);                     // obtain the norm of the gradient of the LevelSet
               VecScale(direccion, vDirection[ip], 1.0 / (gradientnorm + 0.0000000001)); // normalize the gradient of the LevelSet

               // Magnitud of the flux
               number flux = Magnitud; // cantidad de concentración que entrará por las paredes

               VecScale(vValue[ip], direccion, flux);
            }

            return; // salir de de evaluate
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

            if (elem_cut(elem, si, vCornerCoords, time) < 0 ) //  outside (<0) 
            {
               for (size_t ip = 0; ip < nip; ++ip)
               {
                  MathVector<dim> Vel;
                  VecSet(Vel, 0);
                  VecScale(vDarcyVel[ip], Vel, 1.0);
               }
               if (!bDeriv || this->zero_derivative())
                  return;
               this->set_zero(vvvDeriv, nip);
               return;
            }

            std::vector<MathVector<dim>> vDirection(nip);
            (*m_spDirection)(&vDirection[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT); // null para no preguntar por derivada
				const number Magnitud = m_spMagnitud; // 0.0005 // cantidad de concentración que entrará por las paredes

            // const MathVector<dim> *vDirection = m_spDirection->values(s);

            for (size_t ip = 0; ip < nip; ++ip)
            {
               //	Calculate Direccion
               MathVector<dim> direccion;                                              // vector to save direction of the velocity
               number gradientnorm = VecLength(vDirection[ip]);                     // obtain the norm of the gradient of the LevelSet
               VecScale(direccion, vDirection[ip], 1.0 / (gradientnorm + 0.0000000001)); // normalize the gradient of the LevelSet

               // Magnitud of the flux
               number flux = Magnitud; // cantidad de concentración que entrará por las paredes

               VecScale(vDarcyVel[ip], direccion, flux);
            }


            //	check if something to do
            if (!bDeriv || this->zero_derivative())
               return;

            //	clear all derivative values
            this->set_zero(vvvDeriv, nip);

         }

      public:

         ///	set LevelSet gradient import
         void set_gradient_stationary_difussion(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
         {
            m_spDirection = data;
         }

         /// Set domain discretization: SmartPtr<extrapol_type> spExtrapol
         void set_domain_discretizacion(SmartPtr<extrapol_type> spExtrapol)
         {
            m_spExtrapolation = spExtrapol;
         }

			///	set magnitud
			void set_magnitud(number data)
			{
				m_spMagnitud = data;
			}

      protected:

         ///	import for LevelSet gradient
         SmartPtr<CplUserData<MathVector<dim>, dim>> m_spDirection;

         SmartPtr<extrapol_type> m_spExtrapolation;

         number m_spMagnitud;
      };

   } // end namespace NeuroGrowth
} // end namespace ug

#endif /* __H__UG__PLUGINS__LEVEL_SET_DARCY_VELOCITY_LINKER__ */
