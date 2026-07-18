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

#ifndef __H__UG__PLUGINS__LEVEL_SET_CONCENTRATION_VELOCITY_LINKER__
#define __H__UG__PLUGINS__LEVEL_SET_CONCENTRATION_VELOCITY_LINKER__

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
		class LSConcentrationDepentVelocity
			 : public StdDataLinker<LSConcentrationDepentVelocity<TDomain, TAlgebra>, MathVector<TDomain::dim>, TDomain::dim>
		{
			//	domain type
			typedef TDomain domain_type;

			//	algebra type
			typedef TAlgebra algebra_type;

			//	world dimension
			static const int dim = domain_type::dim;

			//	Base class type
			typedef StdDataLinker<LSConcentrationDepentVelocity<domain_type, algebra_type>, MathVector<dim>, dim> base_type;

			//	extrapolation type
			typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

		public:
			LSConcentrationDepentVelocity(SmartPtr<extrapol_type> spExtrapol) : 
				m_spExtrapolation(spExtrapol),
				m_spCalcium(NULL), m_spDCalcium(NULL),
				m_spTubuline(NULL), m_spDTubuline(NULL),
				m_spMAPu(NULL), m_spDMAPu(NULL),
				m_spMAPb(NULL), m_spDMAPb(NULL),
				m_spMAPp(NULL), m_spDMAPp(NULL),
				m_spCurvature(NULL), m_spDCurvature(NULL),
				m_spLevelSetGrad(NULL) , m_spDLevelSetGrad(NULL),
				m_spMinimo(NULL), m_spMaximo(NULL),
				m_spDirection(NULL) , m_spDDirection(NULL),
				m_spminimoTubulin(NULL), m_spminimoInhibition(NULL),
				m_spminimoInhSign(NULL),
				m_spMagnitudVelocity(NULL),
				m_spInhibitionGrad(NULL), m_spDInhibitionGrad(NULL), // poner la gradiente de la inhibicion
				m_spInhibition(NULL), m_spDInhibition(NULL) // poner la inhibicion

			{
				//	this linker needs exactly five input
				this->set_num_input(10);
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
				UG_THROW("LSConcentrationDepentVelocity: Element is necessary for the evaluation.");
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
				if (elem_cut(elem, si, vCornerCoords, time) < 0) // if above the free surface
				{
					for (size_t ip = 0; ip < nip; ++ip)
						vValue[ip] = 0.0;
					return;
				}

				std::vector<number> vCalcium(nip);
				std::vector<number> vTubuline(nip);
				std::vector<number> vMAPu(nip);
				std::vector<number> vMAPb(nip);
				std::vector<number> vMAPp(nip);
				std::vector<number> vInhibition(nip);

				std::vector<number> vCurvature(nip);
           		std::vector<MathVector<dim>> vLevelSetGrad(nip);
				std::vector<MathVector<dim>> vDirection(nip);
				std::vector<MathVector<dim>> vInhibitionGrad(nip);

				// el problema es el mapeo dado que tienen diferente tamaño
				// vCalcium pide 3 -- se verifico que vGlobIP
				// m_spCalcium solo te entrega 2

				(*m_spCalcium)(&vCalcium[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				(*m_spTubuline)(&vTubuline[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				(*m_spMAPu)(&vMAPu[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				(*m_spMAPb)(&vMAPb[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				(*m_spMAPp)(&vMAPp[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            	(*m_spCurvature)(&vCurvature[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            	(*m_spLevelSetGrad)(&vLevelSetGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            	(*m_spDirection)(&vDirection[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				(*m_spInhibition)(&vInhibition[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				(*m_spInhibitionGrad)(&vInhibitionGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);



				const number minimo = m_spMinimo; // 5.5
				const number maximo = m_spMaximo; // 18.5

				const number minimoTubulin = m_spminimoTubulin; //  0.6
				const number minimoInhibition = m_spminimoInhibition; //  0.6
				const number minimoInhSign = m_spminimoInhSign; // 10
				const number MagnitudVelocity = m_spMagnitudVelocity; // 10

				// mapb -> va a ser level-set-function para poder filtrar valores 


				// Marcar los puntos que no son elegibles para tener valor negatico de las intersecciones de las ramas
				std::vector<bool> vMarkedForGrowth(nip, false);


				for (size_t ip = 0; ip < nip; ++ip)
				{
					if (vCurvature[ip] > minimo && vCurvature[ip] < maximo ) // && !vMarkedForGrowth[ip])
               		{
						// normaliyar la gradinete de la concentrqacion de la molecula imaginaria 
						if (vTubuline[ip] > minimoTubulin) //&& vInhibition[ip] <= 0.0008)
						//if (vTubuline[ip] > minimoTubulin)
						{
							//	Calculate Direccion = normal
							MathVector<dim> DirectionGrowth;                                              // vector to save direction of the velocity
							number gradientnorm = VecLength(vDirection[ip]);                     // obtain the norm of the gradient of the LevelSet
							VecScale(DirectionGrowth, vDirection[ip], 1.0 / (gradientnorm + 0.0000001)); // normalize the gradient of the LevelSet

							// Magnitud of the growth
							number growth = MagnitudVelocity * vTubuline[ip]; // cantidad de concentración que entrará por las paredes
							number magnitud = growth * vMAPb[ip]; // cantidad de concentración que entr

							// obtener la velocidad : magnitud . direction
							VecScale(vValue[ip], DirectionGrowth, magnitud);






							if (vInhibition[ip] >= minimoInhibition || vInhibition[ip] >= minimoInhSign) // vInhibition is more than one of the two values
							{
								//  cambiar signo de gradiente de inhibicion: 
								//MathVector<dim> InhibitionInverse;
								//VecScale(InhibitionInverse, vInhibitionGrad[ip], 1.0);
						
								// 1. Calcular la norma del gradiente de inhibición
								number inhibicionNorm = VecLength(vInhibitionGrad[ip]); // Magnitud de grad(inh)

								// 2. Calcular el cuadrado de la norma del gradiente de inhibición
								number inhibicionNorm2 = inhibicionNorm * inhibicionNorm + 0.0000001; // /grad(inh)/^2, evitando divisiones por cero

								// 3. Calcular el producto escalar grad(inh) * V
								number gradInh_V = VecDot(vInhibitionGrad[ip], vValue[ip]); // Producto escalar grad(inh) * V

								// valor absoluto
								gradInh_V = std::abs(gradInh_V);

								// 4. Escalar el gradiente de inhibición con el producto escalar y dividir por /grad(inh)/^2
								MathVector<dim> velInhibicion;
								VecScale(velInhibicion, vInhibitionGrad[ip], gradInh_V / inhibicionNorm2);

								// 5. Calcular la velocidad final: Vf = V - velInhibicion
								MathVector<dim> velFinal;

								//number factor = 1; //0.00013 / vInhibition[ip];
								//VecScale(velInhibicion, velInhibicion, factor);
								//VecScale(velInhibicion, velInhibicion, inhibicionNorm);

								//VecScale(vValue[ip], vValue[ip], 1.0 - factor);
								
								//VecScale(velFinal, velFinal, inhibicionNorm);


								
								/// ELEGIT EL SIGNO (+ or -), esto depende de los 2 umbrales (minimoInhibition, minimoInhSign)
								if (vInhibition[ip] >= minimoInhibition && vInhibition[ip] >= minimoInhSign)
								{
									if (minimoInhibition >= minimoInhSign) {
										VecAdd(velFinal, vValue[ip], velInhibicion);
										VecScale(vValue[ip], velInhibicion, 1.0); // make  growth with correction 
									} else {
										VecSubtract(velFinal, vValue[ip], velInhibicion);
										VecScale(vValue[ip], velFinal, -1.0); // make  retraction
									}
								}
								else if (vInhibition[ip] >= minimoInhibition)
								{
									VecAdd(velFinal, vValue[ip], velInhibicion);
									VecScale(vValue[ip], velInhibicion, 1.0); // make  growth with correction 
								}
								else if (vInhibition[ip] >= minimoInhSign)
								{
									VecSubtract(velFinal, vValue[ip], velInhibicion);
									VecScale(vValue[ip], velFinal, -1.0); // make retraction
								}

							}

						}
						else
						{
							//	Calculate Direccion = normal
							MathVector<dim> DirectionGrowth;                                              // vector to save direction of the velocity
							number gradientnorm = VecLength(vDirection[ip]);                     // obtain the norm of the gradient of the LevelSet
							VecScale(DirectionGrowth, vDirection[ip], 1.0 / (gradientnorm + 0.0000001)); // normalize the gradient of the LevelSet

							//	The sleecte top dont have the enough tubulin, but to mantain their position, we will put a small growth 
							// Magnitud to mantain the position	of the top
							number mag = 0; //0.9; 

							VecScale(vValue[ip], DirectionGrowth, mag);

						}
						// **Marcar este punto como no elegible para branching**
    					vMarkedForGrowth[ip] = false; 
					}
					else if (vCurvature[ip] < -10 && vCurvature[ip] > -35 && VecLength(vValue[ip]) == 0 && vTubuline[ip] < minimoTubulin )
					{
						//	Calculate Direccion
						MathVector<dim> direccion;                                              // vector to save direction of the velocity
						number gradientnorm = VecLength(vLevelSetGrad[ip]);                     // obtain the norm of the gradient of the LevelSet
						VecScale(direccion, vLevelSetGrad[ip], 1.0 / (gradientnorm + 0.00001)); // normalize the gradient of the LevelSet

						// Magnitud of the flux
						number flux = 0; // -0.45 cantidad de concentración que entrará por las paredes

						VecScale(vValue[ip], direccion, flux);

						// **Marcar este punto como no elegible para branching**
    					vMarkedForGrowth[ip] = true; 
					}
					else // This are areas not to be selected for growth or retraction	
					{
						//	Variables
						MathVector<dim> Vel;
						VecSet(Vel, 0);
						VecScale(vValue[ip], Vel, 1.0);

						vMarkedForGrowth[ip] = false;
					}
					

					/// BRANCHING PROCESS : when the MAPp is more 


					
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
				if (elem_cut(elem, si, vCornerCoords, time) < 0) // if above the free surface
				{
					for (size_t ip = 0; ip < nip; ++ip)
						vDarcyVel[ip] = 0.0;
					if (!bDeriv || this->zero_derivative())
						return;
					this->set_zero(vvvDeriv, nip);
					return;
				}

				//	get the data of the ip series
				const number *vCalcium = m_spCalcium->values(s);
				const number *vTubuline = m_spTubuline->values(s);
				const number *vMAPu = m_spMAPu->values(s);
				const number *vMAPb = m_spMAPb->values(s);
				const number *vMAPp = m_spMAPp->values(s);
				const number *vInhibition = m_spInhibition->values(s);

            	const number *vCurvature = m_spCurvature->values(s);
				const MathVector<dim>* vLevelSetGrad = m_spLevelSetGrad->values(s);
				const MathVector<dim>* vDirection = m_spDirection->values(s);
				const MathVector<dim>* vInhibitionGrad = m_spInhibitionGrad->values(s);

				const number minimo = m_spMinimo; // 5.5
				const number maximo = m_spMaximo; // 18.5

				const number minimoTubulin = m_spminimoTubulin; //  0.6

				const number minimoInhibition = m_spminimoInhibition; //  0.6
				const number minimoInhSign = m_spminimoInhSign; // 10

				const number MagnitudVelocity = m_spMagnitudVelocity; //10

				std::vector<bool> vMarkedForGrowth(nip, false); // Inicializa en false para todos los puntos

				for (size_t ip = 0; ip < nip; ++ip)
				{
					if (vCurvature[ip] > minimo && vCurvature[ip] < maximo ) //&& !vMarkedForGrowth[ip])
			   		{
						if (vTubuline[ip] > minimoTubulin) //&& vInhibition[ip] <= 0.0008)
						{
							//	Calculate Direccion = DirectionGrowth
							MathVector<dim> DirectionGrowth;                                              // vector to save direction of the velocity
							number gradientnorm = VecLength(vDirection[ip]);                     // obtain the norm of the gradient of the LevelSet
							VecScale(DirectionGrowth, vDirection[ip], 1.0 / (gradientnorm + 0.0000001)); // normalize the gradient of the LevelSet

							// Magnitud of the growth
							number growth = MagnitudVelocity * vTubuline[ip]; // cantidad de concentración que entrará por las paredes
							number magnitud = growth * vMAPb[ip]; // cantidad de concentración que entr

							// obtener la velocidad : magnitud . direction
							VecScale(vDarcyVel[ip], DirectionGrowth, magnitud);


							if (vInhibition[ip] >= minimoInhibition || vInhibition[ip] >= minimoInhSign) // vInhibition is more than one of the two values
							{
								//  cambiar signo de gradiente de inhibicion: 
								//MathVector<dim> InhibitionInverse;
								//VecScale(InhibitionInverse, vInhibitionGrad[ip], 1.0);
							
								// 1. Calcular la norma del gradiente de inhibición
								number inhibicionNorm = VecLength(vInhibitionGrad[ip]) ; // Magnitud de grad(inh)

								// 2. Calcular el cuadrado de la norma del gradiente de inhibición
								number inhibicionNorm2 = inhibicionNorm * inhibicionNorm + 0.0000001; // /grad(inh)/^2, evitando divisiones por cero

								// 3. Calcular el producto escalar grad(inh) * V
								number gradInh_V = VecDot(vInhibitionGrad[ip], vDarcyVel[ip]); // Producto escalar grad(inh) * V
								// valor absoluto
								gradInh_V = std::abs(gradInh_V);

								// 4. Escalar el gradiente de inhibición con el producto escalar y dividir por /grad(inh)/^2
								MathVector<dim> velInhibicion;
								VecScale(velInhibicion, vInhibitionGrad[ip], gradInh_V / inhibicionNorm2);


								// 5. Calcular la velocidad final: Vf = V - velInhibicion
								MathVector<dim> velFinal;

								//number factor = 1;//0.00013 / vInhibition[ip];
								//VecScale(velInhibicion, velInhibicion, factor);
								//VecScale(velInhibicion, velInhibicion, inhibicionNorm);

								//VecScale(vDarcyVel[ip], vDarcyVel[ip], 1.0 - factor);
								
								//VecScale(velFinal, velFinal, inhibicionNorm);


								/// ELEGIT EL SIGNO (+ or -), esto depende de los 2 umbrales (minimoInhibition, minimoInhSign)
								if (vInhibition[ip] >= minimoInhibition && vInhibition[ip] >= minimoInhSign)
								{
									if (minimoInhibition >= minimoInhSign) {
										VecAdd(velFinal, vDarcyVel[ip], velInhibicion);
										VecScale(vDarcyVel[ip], velInhibicion, 1.0); // make  growth with correction 
									} else {
										VecSubtract(velFinal, vDarcyVel[ip], velInhibicion);
										VecScale(vDarcyVel[ip], velFinal, -1.0); // make  retraction
									}
								}
								else if (vInhibition[ip] >= minimoInhibition)
								{
									VecAdd(velFinal, vDarcyVel[ip], velInhibicion);
									VecScale(vDarcyVel[ip], velInhibicion, 1.0); // make  growth with correction 
								}
								else if (vInhibition[ip] >= minimoInhSign)
								{
									VecSubtract(velFinal, vDarcyVel[ip], velInhibicion);
									VecScale(vDarcyVel[ip], velFinal, -1.0); // make retraction
								}

							}
						
						
						}
						else
						{
							//	Calculate Direccion = normal
							MathVector<dim> DirectionGrowth;                                              // vector to save direction of the velocity
							number gradientnorm = VecLength(vDirection[ip]);                     // obtain the norm of the gradient of the LevelSet
							VecScale(DirectionGrowth, vDirection[ip], 1.0 / (gradientnorm + 0.0000001)); // normalize the gradient of the LevelSet

							//	The sleecte top dont have the enough tubulin, but to mantain their position, we will put a small growth 
							// Magnitud to mantain the position	of the top
							number mag = 0;  //0.9

							VecScale(vDarcyVel[ip], DirectionGrowth, mag);
						}
						// **Marcar este punto como no elegible para branching**
    					vMarkedForGrowth[ip] = false; 
					}
 					else if (vCurvature[ip] < -10 && vCurvature[ip] > -35 && VecLength(vDarcyVel[ip]) == 0  && vTubuline[ip] < minimoTubulin )
					{
						//	Calculate Direccion
						MathVector<dim> direccion;                                              // vector to save direction of the velocity
						number gradientnorm = VecLength(vLevelSetGrad[ip]);                     // obtain the norm of the gradient of the LevelSet
						VecScale(direccion, vLevelSetGrad[ip], 1.0 / (gradientnorm + 0.00001)); // normalize the gradient of the LevelSet

						// Magnitud of the flux
						number flux = -0; //-0.45  cantidad de concentración que entrará por las paredes
						VecScale(vDarcyVel[ip], direccion, flux);
						// **Marcar este punto como no elegible para branching**
    					vMarkedForGrowth[ip] = true; 
					} 
					else
					{
						MathVector<dim> Vel;
						VecSet(Vel, 0);
						VecScale(vDarcyVel[ip], Vel, 1.0);
						vMarkedForGrowth[ip] = false;
					}


				}

				//	check if something to do
					if(!bDeriv || this->zero_derivative()) return;

				//	clear all derivative values
					this->set_zero(vvvDeriv, nip);
			}

		public:
			///	set Calcium import
			void set_calcium(SmartPtr<CplUserData<number, dim>> data)
			{
				m_spCalcium = data;
				m_spDCalcium = data.template cast_dynamic<DependentUserData<number, dim>>();
				base_type::set_input(_RHO_, data, data);
			}
			void set_calcium(number val)
			{
				set_calcium(make_sp(new ConstUserNumber<dim>(val)));
			}

			///	set tubuline import
			void set_tubuline(SmartPtr<CplUserData<number, dim>> data)
			{
				m_spTubuline = data;
				m_spDTubuline = data.template cast_dynamic<DependentUserData<number, dim>>();
				base_type::set_input(_TUB_, data, data);
			}
			void set_tubuline(number val)
			{
				set_tubuline(make_sp(new ConstUserNumber<dim>(val)));
			}

			/// set MAPu import
			void set_MAPu(SmartPtr<CplUserData<number, dim>> data)
			{
				m_spMAPu = data;
				m_spDMAPu = data.template cast_dynamic<DependentUserData<number, dim>>();
				base_type::set_input(_U_, data, data);
			}
			void set_MAPu(number val)
			{
				set_MAPu(make_sp(new ConstUserNumber<dim>(val)));
			}

			/// set MAPb import
			void set_MAPb(SmartPtr<CplUserData<number, dim>> data)
			{
				m_spMAPb = data;
				m_spDMAPb = data.template cast_dynamic<DependentUserData<number, dim>>();
				base_type::set_input(_P_, data, data);
			}
			void set_MAPb(number val)
			{
				set_MAPb(make_sp(new ConstUserNumber<dim>(val)));
			}

			/// set MAPp import
			void set_MAPp(SmartPtr<CplUserData<number, dim>> data)
			{
				m_spMAPp = data;
				m_spDMAPp = data.template cast_dynamic<DependentUserData<number, dim>>();
				base_type::set_input(_B_, data, data);
			}
			void set_MAPp(number val)
			{
				set_MAPp(make_sp(new ConstUserNumber<dim>(val)));
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

			///	set LevelSet gradient import
			void set_LevelSet_gradient(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
			{
				m_spLevelSetGrad = data;
				m_spDLevelSetGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim>>();
				base_type::set_input(_LG_, data, data);
			}


			///	set the direction : this will use the gradient of a imagina
			void set_gradient_stationary_difussion(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
			{
				m_spDirection = data;
				m_spDDirection = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim>>();
				base_type::set_input(_DIR_, data, data);
			}


			///	set Inhibition import
			void set_Inhibition(SmartPtr<CplUserData<number, dim>> data)
			{
				m_spInhibition = data;
				m_spDInhibition = data.template cast_dynamic<DependentUserData<number, dim>>();
				base_type::set_input(_INH_, data, data);
			}

			void set_Inhibition(number val)
			{
				set_Inhibition(make_sp(new ConstUserNumber<dim>(val)));
			}

			///	set Inhibition gradient import
			void set_Inhibition_gradient(SmartPtr<CplUserData<MathVector<dim>, dim>> data)
			{
				m_spInhibitionGrad = data;
				m_spDInhibitionGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim>>();
				base_type::set_input(_INHG_, data, data);
			}


			///	set LevelSet gradient import
			void set_interval_min(number data)
			{
				m_spMinimo = data;
			}

			///	set LevelSet gradient import
			void set_interval_max(number data)
			{
				m_spMaximo = data;
			}

			///	set LevelSet gradient import
			void set_interval_min_Calcium(number data)
			{
				m_spminimoTubulin = data;
			}


			///	set th evalue of minimum the inhibition
			void set_interval_min_Inhibition(number data)
			{
				m_spminimoInhibition = data;
			}

			/// set the minimun value to change the sign of the inhibition
			void set_interval_min_Inhibition_sign(number data)
			{
				m_spminimoInhSign = data;
			}


			///	set LevelSet gradient import
			void set_magnitud_velocity(number data)
			{
				m_spMagnitudVelocity = data;
			}

		protected:
			///	extrapolation by the level-set function
			SmartPtr<extrapol_type> m_spExtrapolation;

			///	import for Calcium
			static const size_t _RHO_ = 0;
			SmartPtr<CplUserData<number, dim>> m_spCalcium;
			SmartPtr<DependentUserData<number, dim>> m_spDCalcium;

			///	import for tubuline
			static const size_t _TUB_ = 1;
			SmartPtr<CplUserData<number, dim>> m_spTubuline;
			SmartPtr<DependentUserData<number, dim>> m_spDTubuline;

			///	import for MAPu
			static const size_t _U_ = 2;
			SmartPtr<CplUserData<number, dim>> m_spMAPu;
			SmartPtr<DependentUserData<number, dim>> m_spDMAPu;

			///	import for MAPb
			static const size_t _B_ = 3;
			SmartPtr<CplUserData<number, dim>> m_spMAPb;
			SmartPtr<DependentUserData<number, dim>> m_spDMAPb;

			///	import for MAPp
			static const size_t _P_ = 4;
			SmartPtr<CplUserData<number, dim>> m_spMAPp;
			SmartPtr<DependentUserData<number, dim>> m_spDMAPp;

			///	import for Curvature
			static const size_t _CUR_ = 5;
			SmartPtr<CplUserData<number, dim>> m_spCurvature;
			SmartPtr<DependentUserData<number, dim>> m_spDCurvature;

			///	import for LevelSet gradient
			static const size_t _LG_ = 6;
			SmartPtr<CplUserData<MathVector<dim>, dim>> m_spLevelSetGrad;
			SmartPtr<DependentUserData<MathVector<dim>, dim>> m_spDLevelSetGrad;

			///	import the gradient of a imaginary molecule = direccion
			static const size_t _DIR_ = 7;
			SmartPtr<CplUserData<MathVector<dim>, dim>> m_spDirection;
			SmartPtr<DependentUserData<MathVector<dim>, dim>> m_spDDirection;

			/// Import for Inhibition
			static const size_t _INH_ = 8;
			SmartPtr<CplUserData<number, dim>> m_spInhibition;
			SmartPtr<DependentUserData<number, dim>> m_spDInhibition;


			///	import for Inhibition gradient (vecotr)
			static const size_t _INHG_ = 9;
			SmartPtr<CplUserData<MathVector<dim>, dim>> m_spInhibitionGrad;
			SmartPtr<DependentUserData<MathVector<dim>, dim>> m_spDInhibitionGrad;

			///	import the values of the interval to select the curvarure
			number m_spMinimo;
			number m_spMaximo;

			///	import the values of minimum value for select the calcium concentration
			number m_spminimoTubulin;

			///	import the values of minimum value for select the inhibition concentration
			number m_spminimoInhibition;

			///	import the values of minimum value for CHANGE THE SIGN OF THE INHIBITION
			number m_spminimoInhSign;

			///	import the values of magnitud
			number m_spMagnitudVelocity;

		};

	} // end namespace NeuroGrowth
} // end namespace ug

#endif /* __H__UG__PLUGINS__LEVEL_SET_DARCY_VELOCITY_LINKER__ */
