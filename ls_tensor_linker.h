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

#ifndef __H__UG__PLUGINS__LEVEL_SET_TENSOR_LINKER__
#define __H__UG__PLUGINS__LEVEL_SET_TENSOR_LINKER__

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
      template <int dim>
      class LSTensorLinker
          : public StdDataLinker<LSTensorLinker<dim>, MathMatrix<dim,dim>, dim>
      {
        //	Base class type
		typedef StdDataLinker< LSTensorLinker<dim>, MathVector<dim>, dim> base_type;


      public:
         LSTensorLinker() : m_spdiffusion_coeff(NULL),
                            m_spGradLevelSet(NULL) 
         {
            //	this linker needs exactly five input
            this->set_num_input(0);
         }

         inline void evaluate(MathMatrix<dim, dim> &value,
                              const MathVector<dim> &globIP,
                              number time, int si) const
         {
            UG_THROW("LSDarcyVelocityLinker: Element is necessary for the evaluation.");
         }



         inline void evaluate(MathMatrix<dim, dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si) const
         {

            MathVector<dim> GradLevelSet;
            (*m_spGradLevelSet)(GradLevelSet, vGlobIP, time, si);

			const number diffusion_coeff = m_spdiffusion_coeff;


            //	Calculate normal
            MathVector<dim> normal;                                              // vector to save direction of the velocity
            number gradientnorm = VecLength(GradLevelSet);                     // obtain the norm of the gradient of the LevelSet
            VecScale(normal, GradLevelSet, 1.0 / (gradientnorm + 0.0000000001)); // normalize the gradient of the LevelSet

            MathMatrix<dim, dim> I;
            MathMatrix<dim, dim> NormalProj;
            MathMatrix<dim, dim> DiffTensor;

            for (size_t i = 0; i < dim; ++i)
            {
                for (size_t j = 0; j < dim; ++j)
                {
                    I[i][j] = (i == j) ? 1.0 : 0.0;
                    NormalProj[i][j] = normal[i] * normal[j];
                    DiffTensor[i][j] = diffusion_coeff * (I[i][j] - NormalProj[i][j]);
                }
            }

            //  save the diffusion tensor
            vValue = DiffTensor;
            
         }


         template <int refDim>
         inline void evaluate(MathMatrix<dim, dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              GridObject *elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<refDim> vLocIP[],
                              const size_t nip,
                              LocalVector *u,
                              const MathMatrix<refDim, dim> *vJT = NULL) const
         {

            std::vector<MathVector<dim>> vGradLevelSet(nip);
            (*m_spGradLevelSet)(&vGradLevelSet[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT);

			const number diffusion_coeff = m_spdiffusion_coeff; // 0.0005 // cantidad de concentración que entrará por las paredes

            for (size_t ip = 0; ip < nip; ++ip)
            {
               //	Calculate normal
               MathVector<dim> normal;                                              // vector to save direction of the velocity
               number gradientnorm = VecLength(vGradLevelSet[ip]);                     // obtain the norm of the gradient of the LevelSet
               VecScale(normal, vGradLevelSet[ip], 1.0 / (gradientnorm + 0.0000000001)); // normalize the gradient of the LevelSet

                MathMatrix<dim, dim> I;
                MathMatrix<dim, dim> NormalProj;
                MathMatrix<dim, dim> DiffTensor;

                for (size_t i = 0; i < dim; ++i)
                {
                    for (size_t j = 0; j < dim; ++j)
                    {
                        I[i][j] = (i == j) ? 1.0 : 0.0;
                        NormalProj[i][j] = normal[i] * normal[j];
                        DiffTensor[i][j] = diffusion_coeff * (I[i][j] - NormalProj[i][j]);
                    }
                }

                //  save the diffusion tensor
                vValue[ip] = DiffTensor;
            }

            return; // salir de de evaluate
         }

         template <int refDim>
         void eval_and_deriv(MathMatrix<dim, dim> vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject *elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector *u,
                             bool bDeriv,
                             int s,
                             std::vector<std::vector<MathMatrix<dim,dim> > > vvvDeriv[],
                             const MathMatrix<refDim, dim> *vJT = NULL) const
         {

        

            std::vector<MathVector<dim>> vGradLevelSet(nip);
            (*m_spGradLevelSet)(&vGradLevelSet[0], vGlobIP, time, si,
                                elem, vCornerCoords, vLocIP, nip, u, vJT); // null para no preguntar por derivada
			
            const number diffusion_coeff = m_spdiffusion_coeff; // 0.0005 // cantidad de concentración que entrará por las paredes

            // const MathVector<dim> *vGradLevelSet = m_spGradLevelSet->values(s);

            for (size_t ip = 0; ip < nip; ++ip)
            {
               //	Calculate normal
               MathVector<dim> normal;                                              // vector to save direction of the velocity
               number gradientnorm = VecLength(vGradLevelSet[ip]);                     // obtain the norm of the gradient of the LevelSet
               VecScale(normal, vGradLevelSet[ip], 1.0 / (gradientnorm + 0.0000000001)); // normalize the gradient of the LevelSet

                MathMatrix<dim, dim> I;
                MathMatrix<dim, dim> NormalProj;
                MathMatrix<dim, dim> DiffTensor;

                for (size_t i = 0; i < dim; ++i)
                {
                    for (size_t j = 0; j < dim; ++j)
                    {
                        I[i][j] = (i == j) ? 1.0 : 0.0;
                        NormalProj[i][j] = normal[i] * normal[j];
                        DiffTensor[i][j] = diffusion_coeff * (I[i][j] - NormalProj[i][j]);
                    }
                }

                //  save the diffusion tensor
                vValue[ip] = DiffTensor;
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
            m_spGradLevelSet = data;
         }


        ///	set diffusion_coeff
        void set_diffusion_coeff(number data)
        {
            m_spdiffusion_coeff = data;
        }

      protected:

         ///	import for LevelSet gradient
         SmartPtr<CplUserData<MathVector<dim>, dim>> m_spGradLevelSet;

         number m_spdiffusion_coeff;
      };

   } // end namespace NeuroGrowth
} // end namespace ug

#endif /* __H__UG__PLUGINS__LEVEL_SET_DARCY_VELOCITY_LINKER__ */
