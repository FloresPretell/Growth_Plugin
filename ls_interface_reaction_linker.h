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

#ifndef __H__UG__PLUGINS__LEVEL_SET_INTERFACE_REACTION_LINKER__
#define __H__UG__PLUGINS__LEVEL_SET_INTERFACE_REACTION_LINKER__

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
      class LSReactionTopBranch
          : public StdDataLinker<LSReactionTopBranch<TDomain, TAlgebra>, number, TDomain::dim>
      {
         //	domain type
         typedef TDomain domain_type;

         //	algebra type
         typedef TAlgebra algebra_type;

         //	world dimension
         static const int dim = domain_type::dim;

         //	Base class type
         typedef StdDataLinker<LSReactionTopBranch<domain_type, algebra_type>,  number, dim> base_type;

         //	extrapolation type
         typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

      public:
         LSReactionTopBranch() : m_spExtrapolation(NULL),
                               m_spFirstTerm(NULL), m_spDFirstTerm(NULL),
                               m_spSecondTerm(NULL), m_spDSecondTerm(NULL),
                               m_spCurvature(NULL), m_spDCurvature(NULL),
                               m_spVelDeformation(NULL), m_spDVelDeformation(NULL),
                               m_spMinimo(NULL), m_spMaximo(NULL),
                               m_spMagnitudFlux(NULL),
                               m_spConstantdiv(NULL), 
                               m_spConstantAddition(NULL),
                               m_spLevelSetGrad(NULL) //, m_spDLevelSetGrad(NULL)
         {
            //	this linker needs exactly five input
            this->set_num_input(4);
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
         inline void evaluate(number vValue,
                              const MathVector<dim> &globIP,
                              number time, int si) const
         {
            UG_THROW("LSDarcyVelocityLinker: Element is necessary for the evaluation.");
         }

         template <int refDim>
         inline void evaluate(number vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              GridObject *elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<refDim> vLocIP[],
                              const size_t nip,
                              LocalVector *u,
                              const MathMatrix<refDim, dim> *vJT = NULL) const
         {

            if (elem_cut(elem, si, vCornerCoords, time) < 0 || elem_cut(elem, si, vCornerCoords, time) > 0) //  outside (<0) and inside (>0) the free surface
            {
               for (size_t ip = 0; ip < nip; ++ip)
               {
                    vValue[ip] = 0.0;
               }
               return; // salir de de evaluate
            }

            std::vector<number> vFirstTerm(nip);
            std::vector<number> vSecondTerm(nip);
            std::vector<number> vCurvature(nip);
            std::vector<number> vVelDeformation(nip);

            std::vector<MathVector<dim>> vLevelSetGrad(nip);

            (*m_spFirstTerm)(&vFirstTerm[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spSecondTerm)(&vSecondTerm[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spCurvature)(&vCurvature[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spVelDeformation)(&vVelDeformation[0], vGlobIP, time, si,
                            elem, vCornerCoords, vLocIP, nip, u, vJT);

            (*m_spLevelSetGrad)(&vLevelSetGrad[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT);
			
            const number MagnitudFlux = m_spMagnitudFlux; // 0.0005 // cantidad de concentración que entrará por las parede
            const number Constantdiv = m_spConstantdiv; // cuadratic term
            const number ConstantAddition = m_spConstantAddition; // adittion term 

            const number minimo = m_spMinimo; // 5.5
            const number maximo = m_spMaximo; // 18.5


            for (size_t ip = 0; ip < nip; ++ip)
            {
               if (vCurvature[ip] > minimo && vCurvature[ip] < maximo)
               {   
                  // Computamos el valor de reacción
                  const number cuadratic = (vFirstTerm[ip] * vFirstTerm[ip]) / 
                                             (Constantdiv + vFirstTerm[ip] * vFirstTerm[ip]);
                  vValue[ip] = vSecondTerm[ip] * MagnitudFlux * cuadratic + ConstantAddition;

                  // Ajustamos el signo dependiendo de vVelDeformation[ip]
                  if (vVelDeformation[ip] < 0.0)
                  {
                        vValue[ip] *= -1;  // Invierte el signo si la velocidad es negativa
                  }
                  else if (vVelDeformation[ip] == 0.0)
                  {
                        vValue[ip] = 0.0;   // Si la velocidad es cero, el valor también lo es
                  }
                  // Si vVelDeformation[ip] > 0, vValue[ip] queda igual (positivo)
               }
               else
               {
                  vValue[ip] = 0.0;  // Si la curvatura está fuera del rango, vValue es 0
               }
            }

            return; // salir de de evaluate
         }

         template <int refDim>
         void eval_and_deriv(number vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject *elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector *u,
                             bool bDeriv,
                             int s,
                             std::vector<std::vector<number> > vvvDeriv[],
                             const MathMatrix<refDim, dim> *vJT = NULL) const
         {

            if (elem_cut(elem, si, vCornerCoords, time) < 0 || elem_cut(elem, si, vCornerCoords, time) > 0) //  outside (<0) and inside (>0) the free surface
            {
               for (size_t ip = 0; ip < nip; ++ip)
               {
                    vValue[ip] = 0.0;
               }
               if (!bDeriv || this->zero_derivative())
                  return;
               this->set_zero(vvvDeriv, nip);
               return;
            }

            //    get the data of the ip series
            const number *vFirstTerm = m_spFirstTerm->values(s);
            const number *vSecondTerm = m_spSecondTerm->values(s);
            const number *vCurvature = m_spCurvature->values(s);
            const number *vVelDeformation = m_spVelDeformation->values(s);


            std::vector<MathVector<dim>> vLevelSetGrad(nip);
            (*m_spLevelSetGrad)(&vLevelSetGrad[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT); // null para no preguntar por derivada
			
            const number MagnitudFlux = m_spMagnitudFlux; // 0.0005 // cantidad de concentración que entrará por las paredes
            const number Constantdiv = m_spConstantdiv; // cuadratic term
            const number ConstantAddition = m_spConstantAddition; // addition term 
           
            const number minimo = m_spMinimo; // 5.5
            const number maximo = m_spMaximo; // 18.5

            // const MathVector<dim> *vLevelSetGrad = m_spLevelSetGrad->values(s);

            for (size_t ip = 0; ip < nip; ++ip)
            {
               if (vCurvature[ip] > minimo && vCurvature[ip] < maximo)
               {   
                  // Computamos el valor de reacción
                  const number cuadratic = (vFirstTerm[ip] * vFirstTerm[ip]) / 
                                             (Constantdiv + vFirstTerm[ip] * vFirstTerm[ip]);
                  vValue[ip] = vSecondTerm[ip] * MagnitudFlux * cuadratic + ConstantAddition;

                  // Ajustamos el signo dependiendo de vVelDeformation[ip]
                  if (vVelDeformation[ip] < 0.0)
                  {
                        vValue[ip] *= -1;  // Invierte el signo si la velocidad es negativa
                  }
                  else if (vVelDeformation[ip] == 0.0)
                  {
                        vValue[ip] = 0.0;   // Si la velocidad es cero, el valor también lo es
                  }
                  // Si vVelDeformation[ip] > 0, vValue[ip] queda igual (positivo)
               }
               else
               {
                  vValue[ip] = 0.0;  // Si la curvatura está fuera del rango, vValue es 0
               }
            }

            //	Compute the derivatives at all ips     //
            /////////////////////////////////////////////

            //	check if something to do
            if (!bDeriv || this->zero_derivative())
               return;

            //	clear all derivative values
            this->set_zero(vvvDeriv, nip);


         }

      public:
         ///	set FirstTerm import
         void set_FirstTerm(SmartPtr<CplUserData<number, dim>> data)
         {
            m_spFirstTerm = data;
            m_spDFirstTerm = data.template cast_dynamic<DependentUserData<number, dim>>();
            base_type::set_input(_RHO_, data, data);
         }

         void set_FirstTerm(number val)
         {
            set_FirstTerm(make_sp(new ConstUserNumber<dim>(val)));
         }

         ///    set SecondTerm import
        void set_SecondTerm(SmartPtr<CplUserData<number, dim>> data)
        {
            m_spSecondTerm = data;
            m_spDSecondTerm = data.template cast_dynamic<DependentUserData<number, dim>>();
            base_type::set_input(_B_, data, data);
        }
        void set_SecondTerm(number val)
        {
            set_SecondTerm(make_sp(new ConstUserNumber<dim>(val)));
        }

         ///	set Curvature import
         void set_Curvature(SmartPtr<CplUserData<number, dim>> data)
         {
            m_spCurvature = data;
            m_spDCurvature = data.template cast_dynamic<DependentUserData<number, dim>>();
            base_type::set_input(_CUR_, data, data);
         }

         void set_Curvature(number val)
         {
            set_Curvature(make_sp(new ConstUserNumber<dim>(val)));
         }


         ///	set Curvature import
         void set_velocity_deformation(SmartPtr<CplUserData<number, dim>> data)
         {
            m_spVelDeformation = data;
            m_spDVelDeformation = data.template cast_dynamic<DependentUserData<number, dim>>();
            base_type::set_input(_VEL_, data, data);
         }

         void set_velocity_deformation(number val)
         {
            set_velocity_deformation(make_sp(new ConstUserNumber<dim>(val)));
         }


         ///	set LevelSet gradient import
         void set_LevelSet_gradient(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
         {
            m_spLevelSetGrad = data;
            // m_spDLevelSetGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim>>();
            // base_type::set_input(_DP_, data, data);
         }

         /// Set domain discretization: SmartPtr<extrapol_type> spExtrapol
         void set_domain_discretizacion(SmartPtr<extrapol_type> spExtrapol)
         {
            m_spExtrapolation = spExtrapol;
         }

         ///	set magnitud del valöor minimo de curvatura import
         void set_interval_min(number data)
         {
            m_spMinimo = data;
         }

         ///	set magnitud del valöor maximo de curvatura import
         void set_interval_max(number data)
         {
            m_spMaximo = data;
         }


        ///	set LevelSet gradient import
        void set_magnitud_influx(number data)
        {
            m_spMagnitudFlux = data;
        }

        void set_constant_div(number data)
        {
            m_spConstantdiv = data;
        }

        void set_constant_addition(number data)
        {
            m_spConstantAddition = data;
        }


      protected:

         ///	import for FirstTerm
         static const size_t _RHO_ = 0;
         SmartPtr<CplUserData<number, dim>> m_spFirstTerm;
         SmartPtr<DependentUserData<number, dim>> m_spDFirstTerm;

        ///	import for SecondTerm
         static const size_t _B_ = 1;
         SmartPtr<CplUserData<number, dim>> m_spSecondTerm;
         SmartPtr<DependentUserData<number, dim>> m_spDSecondTerm;

         ///	import for Curvature
         static const size_t _CUR_ = 2;
         SmartPtr<CplUserData<number, dim>> m_spCurvature;
         SmartPtr<DependentUserData<number, dim>> m_spDCurvature;

         ///	import for Velocity Deformation
         static const size_t _VEL_ = 3;
         SmartPtr<CplUserData<number, dim>> m_spVelDeformation;
         SmartPtr<DependentUserData<number, dim>> m_spDVelDeformation;


         ///	import for LevelSet gradient
         SmartPtr<CplUserData<MathVector<dim>, dim>> m_spLevelSetGrad;

         ///	extrapolation by the level-set function
         SmartPtr<extrapol_type> m_spExtrapolation;

         ///	import the values of the interval to select the curvarure
         number m_spMinimo;
         number m_spMaximo;

         // import magnitud of the reaction (constant)
         number m_spMagnitudFlux;

         // import magnitud of the reaction cuadratic term
         number m_spConstantdiv;

         // import magnitud of the reaction addition term
         number m_spConstantAddition;





      };

   } // end namespace NeuroGrowth
} // end namespace ug

#endif /* __H__UG__PLUGINS__LEVEL_SET_INTERFACE_REACTION_LINKER__ */
