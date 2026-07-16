/*
 * LSSignedNormal — the signed unit normal velocity  w = sign_eps(phi) * grad(phi)/|grad(phi)|.
 *
 * This is the ONE velocity that drives BOTH the reinitialisation and the velocity-extension
 * sweeps of the level-set growth model (the two steps that used the hanging-blind
 * HiResFluxBasedLSM). Feeding this linker as the velocity of the hanging-safe
 * FV1_Convectionhang lets those sweeps run on the adaptive (hanging) mesh:
 *   reinit:     d phi/d tau + w.grad(phi) = sign(phi)     (source = sign, -> |grad phi| -> 1)
 *   extension:  d V  /d tau + w.grad(V)   = 0             (carry V off the interface along normals)
 * The velocity is lagged (recomputed each pseudo-step from the current phi), so the linker
 * only needs its value (zero derivatives), exactly like LSTensorLinker.
 *
 * sign_eps(phi) = phi / sqrt(phi^2 + eps^2) is the smoothed sign (Sussman-Smereka-Osher); eps
 * is a small length ~ h so the interface band is not over-driven. Placed in Nicole's
 * Growth_Plugin/Refinement_LS (may later move to Dmitry's LevelSet — his call).
 *
 * Mirrors the structure of ls_tensor_linker.h (Dmitry Logashenko / Andreas Vogel).
 */
#ifndef __H__UG__PLUGINS__REFINEMENT_LS__LS_SIGNED_NORMAL_LINKER__
#define __H__UG__PLUGINS__REFINEMENT_LS__LS_SIGNED_NORMAL_LINKER__

#include "lib_disc/spatial_disc/user_data/linker/linker.h"
#include "lib_disc/spatial_disc/dom_disc_embb.h"
#include <cmath>

namespace ug
{
   namespace NeuroGrowth
   {
      template <int dim>
      class LSSignedNormal
          : public StdDataLinker<LSSignedNormal<dim>, MathVector<dim>, dim>
      {
         typedef StdDataLinker<LSSignedNormal<dim>, MathVector<dim>, dim> base_type;

      public:
         LSSignedNormal() : m_spGradPhi(NULL), m_spPhi(NULL), m_eps(1e-2)
         {
            this->set_num_input(0);   // imports set via setters (like LSTensorLinker), value-only
         }

      private:
         //	w = sign_eps(phi) * grad / (|grad| + tiny)
         inline void signed_normal(MathVector<dim>& w, const MathVector<dim>& grad, number phi) const
         {
            const number gn = VecLength(grad);
            const number s  = phi / std::sqrt(phi * phi + m_eps * m_eps);   // smoothed sign
            VecScale(w, grad, s / (gn + 1e-12));
         }

      public:
         inline void evaluate(MathVector<dim>& value,
                              const MathVector<dim>& globIP,
                              number time, int si) const
         {
            UG_THROW("LSSignedNormal: Element is necessary for the evaluation.");
         }

         inline void evaluate(MathVector<dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si) const
         {
            MathVector<dim> grad; number phi = 0.0;
            (*m_spGradPhi)(grad, vGlobIP[0], time, si);
            (*m_spPhi)(phi, vGlobIP[0], time, si);
            signed_normal(vValue[0], grad, phi);
         }

         template <int refDim>
         inline void evaluate(MathVector<dim> vValue[],
                              const MathVector<dim> vGlobIP[],
                              number time, int si,
                              GridObject* elem,
                              const MathVector<dim> vCornerCoords[],
                              const MathVector<refDim> vLocIP[],
                              const size_t nip,
                              LocalVector* u,
                              const MathMatrix<refDim, dim>* vJT = NULL) const
         {
            std::vector<MathVector<dim> > vGrad(nip);
            std::vector<number>          vPhi(nip);
            (*m_spGradPhi)(&vGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spPhi)   (&vPhi[0],  vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            for (size_t ip = 0; ip < nip; ++ip)
               signed_normal(vValue[ip], vGrad[ip], vPhi[ip]);
         }

         template <int refDim>
         void eval_and_deriv(MathVector<dim> vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si,
                             GridObject* elem,
                             const MathVector<dim> vCornerCoords[],
                             const MathVector<refDim> vLocIP[],
                             const size_t nip,
                             LocalVector* u,
                             bool bDeriv,
                             int s,
                             std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
                             const MathMatrix<refDim, dim>* vJT = NULL) const
         {
            std::vector<MathVector<dim> > vGrad(nip);
            std::vector<number>          vPhi(nip);
            (*m_spGradPhi)(&vGrad[0], vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            (*m_spPhi)   (&vPhi[0],  vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
            for (size_t ip = 0; ip < nip; ++ip)
               signed_normal(vValue[ip], vGrad[ip], vPhi[ip]);

            //	velocity is lagged (recomputed each pseudo-step) -> no derivative w.r.t. the unknown
            if (!bDeriv || this->zero_derivative()) return;
            this->set_zero(vvvDeriv, nip);
         }

      public:
         ///	set the level-set gradient import (e.g. GridFunctionGradientData(phi,"phi"))
         void set_levelset_gradient(SmartPtr<CplUserData<MathVector<dim>, dim> > data) { m_spGradPhi = data; }
         ///	set the level-set value import (e.g. GridFunctionNumberData(phi,"phi")) — used for the sign
         void set_levelset_value(SmartPtr<CplUserData<number, dim> > data)             { m_spPhi = data; }
         ///	smoothing length of sign_eps(phi) (~ h); larger = gentler near the interface
         void set_eps(number eps) { m_eps = eps; }

      protected:
         SmartPtr<CplUserData<MathVector<dim>, dim> > m_spGradPhi;   // grad(phi)
         SmartPtr<CplUserData<number, dim> >          m_spPhi;       // phi (for the sign)
         number                                       m_eps;        // smoothed-sign length
      };

   } // namespace NeuroGrowth
} // namespace ug

#endif
