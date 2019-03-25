/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Authors: Christian Wehner, Dmitry Logashenko
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

/*
 * High-resolution flux-based level set method: Implementation
 */

// ug4 headers
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_grid/algorithms/attachment_util.h"
#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#endif

namespace ug{
namespace LevelSet{

// Define for the classic ghost fluid method
#define LS_CLASSIC_GHOST_FLUID

/**
 * get the sign of the LSF in an element:
 * 0 if there both the signes,
 * -1 if negative at all the corners
 * 1 if positive at all the corners
 */
template<typename TGridFunction>
int HiResFluxBasedLSM<TGridFunction>::lsf_sign
(
	size_t noc, ///< number of corners
	number lsf [] ///< corner values of the LSF
)
{
	int pos = 0, neg = 0;
	for (size_t co = 0; co < noc; co++)
	{
		if (lsf [co] >= - lsf_threshold () && lsf [co] <= lsf_threshold ())
			return 0; // we consider such elements as intersected; this is important for computation of the gradients
		if (lsf [co] < 0) neg = 1;
		else pos = 1;
	}
	return pos - neg;
}

/**
 * Marks corners of intersected elements
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::mark_CoIE
(
	grid_type& grid, ///< the grid
	ABool& aCoIE, ///< the attachment for the marks
	t_aaCoIE& aaCoIE ///< accessor for the attachment
)
{
//	set the default value (false everywhere)
	SetAttachmentValues (aaCoIE, grid.vertices_begin (), grid.vertices_end (), false);
	
//	if no LSF then nothing to do
	if (! m_spLSF.valid ()) return;
	
	LocalIndices locInd;
	LocalVector locLSF;
	
//	loop the vertices
	ElemIterator iterEnd = m_spLSF->template end<ElemType> ();
	for (ElemIterator iter = m_spLSF->template begin<ElemType> (); iter != iterEnd; ++iter)
	{
	//	local values of the LSF
		number lsf [maxNumCo];
		ElemType* elem = *iter;
		size_t noc = elem->num_vertices ();
		m_spLSF->indices (elem, locInd);
		locLSF.resize (locInd);
		GetLocalVector (locLSF, *m_spLSF);
		for (size_t i = 0; i < noc; i++)
			lsf[i] = locLSF (0, i);
	//	check the sign of the lsf
		if (lsf_sign (noc, lsf) == 0)
			for (size_t i = 0; i < elem->num_vertices (); i++)
				aaCoIE[elem->vertex (i)] = true;
	}
	
#	ifdef UG_PARALLEL
	AttachmentAllReduce<Vertex> (grid, aCoIE, PCL_RO_LOR);
#	endif
}

/**
 * extrapolate a grid function according to the lsf
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::extrapolate_by_lsf
(
	const CplUserData<number,dim> * if_val_data, ///< computes the values at the interface (if not NOLL)
	int si, ///< subset index (only used if if_val_data != NULL)
	DimFV1Geometry<dim> & geo, ///< FV geometry
	number sol [], ///< corner values of the function to extrapolate
	number lsf [], ///< corner values of the LSF
	size_t base, ///< the base corner
	number ext [] ///< extrapolated values (at all the corners)
)
{
	const MathVector<dim>* co_coord = geo.corners ();
	
	for (size_t co = 0; co < geo.num_scv (); co++)
		if (lsf [co] * lsf [base] > 0)
			ext[co] = sol[co]; /* take the original value */
		else /* extrapolate */
		{
			number interface_val;
			if (if_val_data == NULL)
				interface_val = 0;
			else
			{
				const MathVector<dim>* co_local = geo.scv_local_ips ();
				MathVector<dim> if_pnt_coord, if_pnt_local;
				number s = lsf[co] / (lsf[co] - lsf[base]);
				VecScaleAdd (if_pnt_coord, s, co_coord[base], 1 - s, co_coord[co]);
				VecScaleAdd (if_pnt_local, s, co_local[base], 1 - s, co_local[co]);
				(* if_val_data) (&interface_val, &if_pnt_coord, m_time, si,
					geo.elem (), co_coord, &if_pnt_local, 1, NULL);
			}
			
			number t = lsf[co] / lsf[base];
			ext[co] = sol[base] * t + interface_val * (1 - t);
		}
}

/**
 * compute the update of the solution through a SCVF and the correction due to the source
 */
template<typename TGridFunction>
inline void HiResFluxBasedLSM<TGridFunction>::sol_update
(
	bool redOrder, ///< whether to reduce the order of the discretization
	const MathVector<dim>& ip, ///< coordinates of the integration point
//	Data at the upwind corner
	const MathVector<dim>& x_up, ///< coordinates of the upwind corner
	number u_up, ///< solution at the upwind corner
	const MathVector<dim>& grad_up, ///< gradient at the upwind corner
	const MathVector<dim>& vel_up, ///< velocity at the upwind corner
//	Data at the downwind corner
	const MathVector<dim>& x_down, ///< coordinates of the upwind corner
	number u_down, ///< solution at the downwind corner
	const MathVector<dim>& grad_down, ///< gradient at the downwind corner
	const MathVector<dim>& vel_down, ///< velocity at the downwind corner
//	Computed update
	number& corr_up, ///< update for the upwind corner
	number& corr_down, ///< update for the downwind corner
	number& src_up, ///< source factor for the upwind corner
	number& src_down ///< source factor at the downwind corner
)
{
	if (! m_firstOrder) // the High-Resolution Flux-Based Level-Set method
	{
	//	due to the convection:
	//	conv_corr = u_{ip}^{n+0.5}, where u_{ip}^{n+0.5} is interpolated along the characteristic
		MathVector<dim> distVec;
		VecSubtract (distVec, ip, x_up);
		corr_up = u_up - 0.5 * m_dt * (grad_up * vel_up);
		if (! redOrder)
			corr_up += grad_up * distVec;
		corr_down = corr_up;
	
	//	due to the divergence:
		if (! m_divFree)
		{
		//	div_corr_up = u_{co_up}^{n+0.5}, where u_{co_up}^{n+0.5} is interpolated along the characteristic
			corr_up -= u_up - 0.5 * m_dt * (grad_up * vel_up);
	
		//	div_corr_down = u_{co_down}^{n+0.5}, where u_{co_up}^{n+0.5} is interpolated along the characteristic
			corr_down -= u_down - 0.5 * m_dt * (grad_down * vel_down);
		}
	}
	else // the classic (first-order) upwind method
	{
	//	due to the convection
		corr_down = corr_up = u_up;
	
	//	due to the divergence
		if (! m_divFree)
		{
			corr_up -= u_up;
			corr_down -= u_down;
		}
	}
	
	if (! m_antiderivSrc) // the usual discretization
		src_up = src_down = 0; // we do not assemble the source factors here
	else
	{ // approximate the mesh size
		MathVector<dim> hVec;
		VecSubtract (hVec, x_down, x_up);
		src_up = 0;
		src_down = vel_down * hVec;
	}
}

/**
 * compute the update of the solution through a BF and the correction due to the source
 */
template<typename TGridFunction>
inline void HiResFluxBasedLSM<TGridFunction>::bnd_sol_update
(
	bool redOrder, ///< whether to reduce the order of the discretization
	const MathVector<dim>& bip, ///< coordinates of the boundary integration point
//	Data at the upwind corner
	const MathVector<dim>& x, ///< coordinates of the boundary corner
	number u, ///< solution at the corner
	const MathVector<dim>& grad, ///< gradient at the corner
	const MathVector<dim>& vel, ///< velocity at the corner
//	Computed update
	number& corr ///< update for the corner
)
{
	if (! m_firstOrder) // the High-Resolution Flux-Based Level-Set method
	{
	//	due to the convection:
	//	conv_corr = u_{bip}^{n+0.5}, where u_{bip}^{n+0.5} is interpolated along the characteristic
		MathVector<dim> distVec;
		VecSubtract (distVec, bip, x);
		corr = u - 0.5 * m_dt * (grad * vel);
		if (! redOrder)
			corr += grad * distVec;
	
	//	due to the divergence
		if (! m_divFree)
		//	div_corr = u_{co}^{n+0.5}, where u_{co}^{n+0.5} is interpolated along the characteristic
			corr -= u - 0.5 * m_dt * (grad * vel);
	}
	else // the classic (first-order) upwind method
	{
	//	due to the convection
		corr = u;
	
	//	due to the divergence
		if (! m_divFree)
			corr -= u;
	}
}

/**
 * get the corner velocity
 * We assume that the element is not intersected by the interface.
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::get_nodal_vel
(
	ElemType* elem, ///< the element to compute the contribution for
	MathVector<dim> co_coord[], ///< coordinates of the corners
	DimFV1Geometry<dim>& geo, ///< structure for the FV geometry
	LocalVector& u, ///< solution on the old time level
	MathVector<dim> grad[], ///< computed gradient at corners
	MathVector<dim> co_vel[], ///< array for the nodal velocity
	int lsf_sign ///< sign (1 or -1) of the lsf in the element (or 0 if no lsf specified!)
)
{
	size_t noc = geo.num_scv ();
	
	const int si = 0; //TODO this should be corrected
	
//	Compute the corner velocity
	if (lsf_sign == 0 && m_imVelocity.valid ()) // we consider an arbitrary vel. field only if there is no LSF
	{
		(*m_imVelocity) (co_vel, geo.scv_global_ips (), m_time, si,
			elem, co_coord, geo.scv_local_ips (), noc, &u);
        if (m_gamma != 1)
        	for (size_t i = 0; i < noc; i++)
        		co_vel[i] *= m_gamma;
	}
	else if (lsf_sign == 0 && m_imNormalVel.valid ()) // normal velocity (only if there is no LSF)
	{
		number normal_vel [maxNumCo];
		
		(*m_imNormalVel) (normal_vel, geo.scv_global_ips (), m_time, si,
			elem, co_coord, geo.scv_local_ips (), noc, &u);
			
    	for (size_t i = 0; i < noc; i++)
    	{
    	    number vnorm = VecLength (grad[i]);
    	    if (vnorm > 1e-15) //TODO: Eliminate the explicit constant
    	    {
    	    	if (lsf_sign < 0)
    	    		vnorm = - vnorm;
    	    	VecScale (co_vel[i], grad[i], m_gamma * normal_vel[i] / vnorm);
    	    }
    	    else
    	    	co_vel[i] = 0;
    	}
	}
	else if (m_delta != 0) // constantly scaled normalized gradient as the velocity
    	for (size_t i = 0; i < noc; i++)
    	{
    	    number vnorm = VecLength (grad[i]);
    	    if (vnorm > 1e-15) //TODO: Eliminate the explicit constant
    	    {
    	    	if (lsf_sign < 0)
    	    		vnorm = - vnorm;
    	    	VecScale (co_vel[i], grad[i], m_delta / vnorm);
    	    }
    	    else
    	    	co_vel[i] = 0;
    	}
    else
		for (size_t i = 0; i < noc; i++) co_vel[i] = 0;
}

/**
 * assemble a usual grid element using upwind
 *
 * This function computes the contribution of the explicit local discretization
 * in one grid element.
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::assemble_element
(
	ElemType* elem, ///< the element to compute the contribution for
	DimFV1Geometry<dim>& geo, ///< structure for the FV geometry
	domain_type& domain, ///< domain
	LocalVector& uOld, ///< solution on the old time level
	t_aaGrad& aaGradient, ///< computed gradient at vertices
	t_aaGrad& aaVelGrad, ///< computed gradient at vertices for the computation of the velocity
	t_aaVol& aaVolume, ///< volumes of the SCVs (assigned to vertices)
	int sign, ///< sign of the LSF in the element (or 0 if no LSF specified)
	t_aaUpd& aaUpdate, ///< to accumulate the advective update
	t_aaUpd* aaSrc ///< to accumulate the factor of the update due to the source (if needed)
)
{
//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor ();
	
//	get vertices and extract corner coordinates
	MathVector<dim> coCoord[maxNumCo];
	Vertex* vVrt[maxNumCo];
	for (size_t i = 0; i < elem->num_vertices (); ++i)
	{
		vVrt[i] = elem->vertex (i);
		coCoord[i] = aaPos[vVrt[i]];
	}

//	update fv geometry
	geo.update (elem, coCoord, domain.subset_handler().get ());
	size_t noc = geo.num_scv ();

//  fill node values and gradients
	number uValue[maxNumCo];
    MathVector<dim> grad[maxNumCo];
    MathVector<dim> vel_grad[maxNumCo];
	for (size_t i=0; i < noc; i++)
	{
		uValue[i] = uOld (0, i);
		grad[i] = aaGradient[vVrt[i]];
		vel_grad[i] = aaVelGrad[vVrt[i]];
	}
	
//	get corner velocity and source
	MathVector<dim> coVelocity[maxNumCo];
	get_nodal_vel (elem, coCoord, geo, uOld, vel_grad, coVelocity, sign);

//	outflow boundary
	bool outBndCo[maxNumCo];
	for (size_t i = 0; i < noc; i++) outBndCo[i] = false;
	for (size_t k = 0; k < m_neumann_sg.size (); k++)
	{
		int si = m_neumann_sg [k];
		for (size_t i = 0; i < geo.num_bf (si); i++)
		{
		// 	get current BF
			const typename DimFV1Geometry<dim>::BF& bf = geo.bf (si, i);
			const size_t nodeID = bf.node_id ();
			
		//	mark the corner
			outBndCo[nodeID] = true;
			
		//	compute values at the bip
			MathVector<dim> bipVelocity;
			bipVelocity = 0;
			for (size_t co = 0; co < noc; co++)
				VecScaleAppend (bipVelocity, bf.shape (co), coVelocity[co]);
			number bipNormalVel = bipVelocity * bf.normal ();
			
		//	assemble the fluxes
			number corr;
			bnd_sol_update (true, bf.global_ip (), coCoord[nodeID], uValue[nodeID],
				grad[nodeID], coVelocity[nodeID], corr);
			aaUpdate[ vVrt[nodeID] ] -= bipNormalVel * corr / aaVolume[ vVrt[nodeID] ];
		
		//	the local Courant-number
			if (m_spCourant.valid ())
			{
				number localCFL = m_dt * bipNormalVel / aaVolume[ vVrt[nodeID] ];
				const LocalIndices& ind = uOld.get_indices ();
				const size_t index = ind.index (0, nodeID), comp = ind.comp (0, nodeID);
				number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
				if (m_divFree)
				{ // we consider the outflow faces
					if (bipNormalVel > 0)
						CFL_entry += localCFL;
				}
				else
				{ // actually, CFL should NOT be regarded here...
					CFL_entry += std::fabs (localCFL);
				}
			}
		}
	}
	
//	fluxes through the inner scvfaces
	for (size_t ip = 0; ip < geo.num_scvf (); ++ip)
	{
	//	get current SCVF
	    const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf (ip);
	    size_t from = scvf.from ();
	    size_t to   = scvf.to ();
	    const MathVector<dim>& ipCoord = scvf.global_ip ();
	    
	//  compute the ip velocity from the corner velocity by the linear interpolation
		MathVector<dim> ipVelocity;
		ipVelocity = 0;
		for (size_t co = 0; co < noc; co++)
			VecScaleAppend (ipVelocity, scvf.shape (co), coVelocity[co]);
	
	//	normal ip-velocity
		number ipNormalVel = ipVelocity * scvf.normal ();
		
	//	upwinding
		size_t up_co, down_co;
	    if (ipNormalVel > 0)
	    {
		    up_co = from; down_co = to;
		}
		else
		{
		    up_co = to; down_co = from; ipNormalVel = - ipNormalVel;
		}
		
	//	assemble the fluxes
		number corr_up, corr_down;
		number src_up, src_down;
		sol_update (outBndCo[up_co], // reduce the order at the outflow boundary
			ipCoord, coCoord[up_co], uValue[up_co], grad[up_co], coVelocity[up_co],
				coCoord[down_co], uValue[down_co], grad[down_co], coVelocity[down_co],
					corr_up, corr_down, src_up, src_down);
		
		aaUpdate[ vVrt[up_co] ] -= ipNormalVel * corr_up / aaVolume[ vVrt[up_co] ];
		aaUpdate[ vVrt[down_co] ] += ipNormalVel * corr_down / aaVolume[ vVrt[down_co] ];
		
		if (aaSrc != NULL)
		{
			(*aaSrc)[ vVrt[up_co] ] -= ipNormalVel * src_up / aaVolume[ vVrt[up_co] ];
			(*aaSrc)[ vVrt[down_co] ] += ipNormalVel * src_down / aaVolume[ vVrt[down_co] ];
		}
		
	// the local Courant-number
		if (m_spCourant.valid ())
		{
			const LocalIndices& ind = uOld.get_indices ();
			size_t index, comp;
			
			if (m_divFree)
			{	// consider the outflow faces:
				number localCFL = m_dt * ipNormalVel / aaVolume[ vVrt[up_co] ];
				index = ind.index (0, up_co); comp = ind.comp (0, up_co);
				number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
				CFL_entry += localCFL;
			}
			else
			{	// consider inflow faces
				number localCFL = m_dt * ipNormalVel / aaVolume[ vVrt[down_co] ];
				index = ind.index (0, down_co); comp = ind.comp (0, down_co);
				number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
				CFL_entry += localCFL;
			}
		}
	}
}

/**
 * get the velocity for a given SCVF in an element intersected by the interface
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::get_scvf_vel_on_if
(
	DimFV1Geometry<dim>& geo, ///< structure for the FV geometry
	const typename DimFV1Geometry<dim>::SCVF& scvf, ///< the SCVF to compute the values for
	number u[], ///< corner values of the potential of the velocity
	MathVector<dim> grad[], ///< computed gradient at corners
	number lsf[], ///< corner values of the LSF
//---- Data computed for the from-corner
	MathVector<dim>& from_co_vel, ///< corner velocity for the from-corner
	number& from_flux, ///< normal velocity at the SCVF for the from-corner
//---- Data computed for the to-corner
	MathVector<dim>& to_co_vel, ///< corner velocity for the to-corner
	number& to_flux ///< normal velocity at the SCVF for the to-corner
)
{
	number ext_u [maxNumCo];
	MathVector<dim> grad_ip;
	number norm, delta;
	size_t co, noc = geo.num_scv ();
	
//	From-corner
	co = scvf.from ();
	delta = (lsf[co] >= 0)? 1 : -1;
	extrapolate_by_lsf (NULL, 0, geo, u, lsf, co, ext_u);
	grad_ip = 0.0;
	for (size_t sh = 0; sh < noc; sh++)
		VecScaleAppend (grad_ip, ext_u[sh], scvf.global_grad (sh));
	norm = VecLength (grad_ip);
	if (norm > 1e-15) //TODO: Eliminate the explicit constant
		from_flux = delta * (grad_ip * scvf.normal ()) / norm;
	else
		from_flux = 0;
	norm = VecLength (grad[co]);
	if (norm > 1e-15) //TODO: Eliminate the explicit constant
		VecScale (from_co_vel, grad[co], delta / norm);
	else
		from_co_vel = 0;
	
//	To-corner
	co = scvf.to ();
	delta = (lsf[co] >= 0)? 1 : -1;
	extrapolate_by_lsf (NULL, 0, geo, u, lsf, co, ext_u);
	grad_ip = 0.0;
	for (size_t sh = 0; sh < noc; sh++)
		VecScaleAppend (grad_ip, ext_u[sh], scvf.global_grad (sh));
	norm = VecLength (grad_ip);
	if (norm > 1e-15) //TODO: Eliminate the explicit constant
		to_flux = delta * (grad_ip * scvf.normal ()) / norm;
	else
		to_flux = 0;
	norm = VecLength (grad[co]);
	if (norm > 1e-15) //TODO: Eliminate the explicit constant
		VecScale (to_co_vel, grad[co], delta / norm);
	else
		to_co_vel = 0;
}

/**
 * get the velocity for a given BF in an element intersected by the interface
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::get_bf_vel_on_if
(
	DimFV1Geometry<dim>& geo, ///< structure for the FV geometry
	const typename DimFV1Geometry<dim>::BF& bf, ///< the BF to compute the values for
	number u[], ///< corner values of the potential of the velocity
	MathVector<dim> grad[], ///< computed gradient at corners
	number lsf[], ///< corner values of the LSF
	MathVector<dim>& co_vel, ///< corner velocity
	number& flux ///< normal velocity at the BF
)
{
	number ext_u [maxNumCo];
	MathVector<dim> grad_ip;
	number norm, delta;
	size_t co, noc = geo.num_scv ();
	
	co = bf.node_id ();
	delta = (lsf[co] >= 0)? 1 : -1;
	extrapolate_by_lsf (NULL, 0, geo, u, lsf, co, ext_u);
	grad_ip = 0.0;
	for (size_t sh = 0; sh < noc; sh++)
		VecScaleAppend (grad_ip, ext_u[sh], bf.global_grad (sh));
	norm = VecLength (grad_ip);
	if (norm > 1e-15) //TODO: Eliminate the explicit constant
		flux = delta * (grad_ip * bf.normal ()) / norm;
	else
		flux = 0;
	norm = VecLength (grad[co]);
	if (norm > 1e-15) //TODO: Eliminate the explicit constant
		VecScale (co_vel, grad[co], delta / norm);
	else
		co_vel = 0;
}

/**
 * assemble an element intersected by the interface
 *
 * This function computes the contribution of the explicit local discretization
 * in one grid element intersected by the interface. It checks whether the element
 * is really intersected. If yes, it computes the contributions and returns 0.
 * If no, it does nothing and returns the sign of the lsf in the element.
 */
template<typename TGridFunction>
int HiResFluxBasedLSM<TGridFunction>::assemble_cut_element
(
	ElemType* elem, ///< the element to compute the contribution for
	DimFV1Geometry<dim>& geo, ///< structure for the FV geometry
	domain_type& domain, ///< domain
	LocalVector& uOld, ///< solution on the old time level
	LocalVector& locLSF, ///< the level-set function that defines the interface
	LocalVector& locVelPot, ///< the "velocity potential"
	t_aaGrad& aaGradient, ///< computed gradient at vertices
	t_aaGrad& aaVelGrad, ///< computed gradient at vertices for the computation of the velocity
	t_aaVol& aaVolume, ///< volumes of the SCVs (assigned to vertices)
	t_aaUpd& aaUpdate, ///< to accumulate the convective update
	t_aaUpd* aaSrc, ///< to accumulate the factors of the update due to the source
	CplUserData<number,dim> * if_val_data, ///< computes the values at the interface (if not NOLL)
	int si ///< subset index (only used if if_val_data != NULL)
)
{
//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor ();
	
//	get the LSF and compute its sign
	number lsf[maxNumCo];
	GetLocalVector (locLSF, *m_spLSF);
	for (size_t i = 0; i < elem->num_vertices (); i++)
		lsf[i] = locLSF (0, i);
	int sign = lsf_sign (elem->num_vertices (), lsf);
	if (sign != 0)
		return sign;

//	get vertices and extract corner coordinates
	MathVector<dim> coCoord[maxNumCo];
	Vertex* vVrt[maxNumCo];
	for (size_t i = 0; i < elem->num_vertices (); i++)
	{
		vVrt[i] = elem->vertex (i);
		coCoord[i] = aaPos[vVrt[i]];
	}
	
//	update fv geometry
	geo.update (elem, coCoord, domain.subset_handler().get ());
	size_t noc = geo.num_scv ();

//  fill nodal values and gradients
	number uValue[maxNumCo];
	number vel_pot[maxNumCo];
    MathVector<dim> grad[maxNumCo];
    MathVector<dim> vel_grad[maxNumCo];
	GetLocalVector (locVelPot, *m_spVelPot);
	for (size_t i = 0; i < noc; i++)
	{
		uValue[i] = uOld (0, i);
		vel_pot[i] = locVelPot (0, i);
		grad[i] = aaGradient[vVrt[i]];
		vel_grad[i] = aaVelGrad[vVrt[i]];
	}
	
//	outflow boundary
	for (size_t k = 0; k < m_neumann_sg.size (); k++)
	{
		int si = m_neumann_sg [k];
		for(size_t i = 0; i < geo.num_bf (si); i++)
		{
		// 	get the BF
			const typename DimFV1Geometry<dim>::BF& bf = geo.bf (si, i);
			const size_t nodeID = bf.node_id ();
			
		//	get the velocity
			MathVector<dim> co_vel;
			number flux;
			get_bf_vel_on_if (geo, bf, vel_pot, vel_grad, lsf, co_vel, flux);
			
		//	assemble the fluxes
			if (m_divFree) // otherwise we do not consider any flux for the first order approx.
				aaUpdate[ vVrt[nodeID] ] -= flux * uValue[nodeID] / aaVolume[ vVrt[nodeID] ]; // first order approximation
		
		// the local Courant-number
			if (m_spCourant.valid ())
			{
				if (m_divFree) // otherwise there are no contribution
				{
					if (flux > 0) // we consider only outflow faces
					{
						number localCFL = m_dt * flux / aaVolume[ vVrt[nodeID] ];
						const LocalIndices& ind = uOld.get_indices ();
						const size_t index = ind.index (0, nodeID), comp = ind.comp (0, nodeID);
						number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
						CFL_entry += localCFL;
					}
				}
			}
		}
	}
	
#	ifndef LS_CLASSIC_GHOST_FLUID

//	fluxes through the inner scv faces
	for (size_t ip = 0; ip < geo.num_scvf (); ip++)
	{
		number corr, src, t;
		
	//	get current SCVF
	    const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf (ip);
	    const MathVector<dim>& ipCoord = scvf.global_ip ();
	    size_t from = scvf.from ();
	    size_t to   = scvf.to ();
	    
	//  compute the velocities
		MathVector<dim> from_co_vel, to_co_vel;
		number from_flux, to_flux;
		get_scvf_vel_on_if (geo, scvf, vel_pot, vel_grad, lsf,
			from_co_vel, from_flux, to_co_vel, to_flux);
	
	//	assemble the flux for the from-corner
		sol_update (false, ipCoord, coCoord[from], uValue[from], grad[from], from_co_vel,
			coCoord[to], uValue[to], grad[to], to_co_vel, corr, t, src, t);
		aaUpdate[ vVrt[from] ] -= from_flux * corr / aaVolume[ vVrt[from] ];
		if (aaSrc != NULL)
			(*aaSrc) [ vVrt[from] ] -= from_flux * src / aaVolume[ vVrt[from] ];
		
	//	assemble the flux for the to-corner
		sol_update (false, ipCoord, coCoord[to], uValue[to], grad[to], to_co_vel,
			coCoord[from], uValue[from], grad[from], from_co_vel, corr, t, src, t);
		aaUpdate[ vVrt[to] ] += to_flux * corr / aaVolume[ vVrt[to] ];
		if (aaSrc != NULL)
			(*aaSrc) [ vVrt[to] ] += to_flux * src / aaVolume[ vVrt[to] ];
		
	// the local Courant-number
		if (m_spCourant.valid ())
		{
			const LocalIndices& ind = uOld.get_indices ();
			size_t index, comp;
			
			if (m_divFree)
			{	// consider the outflow faces:
				number localCFL = m_dt * std::abs (from_flux) / aaVolume[ vVrt[from] ];
				index = ind.index (0, from); comp = ind.comp (0, from);
				number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
				CFL_entry += localCFL;
			}
			else
			{	// consider inflow faces
				number localCFL = m_dt * std::abs (to_flux) / aaVolume[ vVrt[to] ];
				index = ind.index (0, to); comp = ind.comp (0, to);
				number& CFL_entry_2 = BlockRef ((* m_spCourant) [index], comp);
				CFL_entry += localCFL;
			}
		}
	}
	
#	else // LS_CLASSIC_GHOST_FLUID

//	loop the base corner for the extrapolation
	for (size_t base_co = 0; base_co < noc; base_co++)
	{
	//	extrapolate the solution
		number extValue [maxNumCo];
		extrapolate_by_lsf (if_val_data, si, geo, uValue, lsf, base_co, extValue);
		
	//	fluxes through the inner scvfaces
		for (size_t ip = 0; ip < geo.num_scvf (); ++ip)
		{
		//	get current SCVF
			const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf (ip);
			size_t from = scvf.from ();
			size_t to   = scvf.to ();
			if (from != base_co && to != base_co)
				continue; // not the current case
			
			const MathVector<dim>& ipCoord = scvf.global_ip ();
		
		//  compute the velocities
			MathVector<dim> from_co_vel, to_co_vel;
			number from_ipNormalVel, to_ipNormalVel;
			get_scvf_vel_on_if (geo, scvf, vel_pot, vel_grad, lsf,
				from_co_vel, from_ipNormalVel, to_co_vel, to_ipNormalVel);
	
		//	normal ip-velocity (depending on the corner we assemble):
			number ipNormalVel = (base_co == from)? from_ipNormalVel : to_ipNormalVel;
		
		//	upwinding
			size_t up_co, down_co;
			MathVector<dim> up_co_vel, down_co_vel;
			if (ipNormalVel > 0)
			{
				up_co = from; down_co = to;
				up_co_vel = from_co_vel; down_co_vel = to_co_vel;
			}
			else
			{
				up_co = to; down_co = from; ipNormalVel = - ipNormalVel;
				up_co_vel = to_co_vel; down_co_vel = from_co_vel;
			}
		
		//	assemble the fluxes
			number corr_up, corr_down;
			number src_up, src_down;
			sol_update (false,
				ipCoord, coCoord[up_co], extValue[up_co], grad[up_co], up_co_vel,
					coCoord[down_co], extValue[down_co], grad[down_co], down_co_vel,
						corr_up, corr_down, src_up, src_down);
			if (base_co == up_co)
			{
				aaUpdate[ vVrt[base_co] ] -= ipNormalVel * corr_up / aaVolume[ vVrt[up_co] ];
				if (aaSrc != NULL)
					(*aaSrc) [ vVrt[base_co] ] -= ipNormalVel * src_up / aaVolume[ vVrt[up_co] ];
			}
			else
			{
				aaUpdate[ vVrt[base_co] ] += ipNormalVel * corr_down / aaVolume[ vVrt[down_co] ];
				if (aaSrc != NULL)
					(*aaSrc) [ vVrt[base_co] ] += ipNormalVel * src_down / aaVolume[ vVrt[down_co] ];
			}
		
		// the local Courant-number
			if (m_spCourant.valid ())
			{
				const LocalIndices& ind = uOld.get_indices ();
				size_t index, comp;
			
				if (m_divFree)
				{	// consider the outflow faces:
					number localCFL = m_dt * ipNormalVel / aaVolume[ vVrt[from] ];
					index = ind.index (0, from); comp = ind.comp (0, from);
					number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
					CFL_entry += localCFL;
				}
				else
				{	// consider inflow faces
					number localCFL = m_dt * ipNormalVel / aaVolume[ vVrt[to] ];
					index = ind.index (0, to); comp = ind.comp (0, to);
					number& CFL_entry = BlockRef ((* m_spCourant) [index], comp);
					CFL_entry += localCFL;
				}
			}
		}
	}
	
#	endif
	
	return 0;
}

/**
 * computes volume of control volumes and saves the result in the grid attachment
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::compute_volumes
(
	TGridFunction& u, ///< grid function to get the grid and the surface view
	DimFV1Geometry<dim>& geo, ///< FV geometry object
	ANumber& aVolume, ///< where to save the volumes
	t_aaVol& aaVolume ///< accessor for the attachment
)
{
//	get domain
	domain_type& domain = *u.domain().get ();

//	get grid of domain
	grid_type& grid = *domain.grid ();

//	initialize attachment value
	SetAttachmentValues (aaVolume, grid.vertices_begin (), grid.vertices_end (), 0);

//	local values
	MathVector<dim> coCoord[maxNumCo];
	Vertex* vVrt[maxNumCo];

//	sum up all contributions of the sub control volumes to one vertex in an attachment
	for (int si = 0; si < domain.subset_handler()->num_subsets (); si++)
	{
	//	loop grid elements of the full dimensionality
		ElemIterator iterEnd = u.template end<ElemType> (si);
		for (ElemIterator iter = u.template begin<ElemType> (si); iter != iterEnd; ++iter)
		{
			ElemType* elem = *iter;

		//	get position accessor
			const position_accessor_type& aaPos = domain.position_accessor ();

		//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices ();
			for (size_t i = 0; i < numVertices; i++)
			{
				vVrt[i] = elem->vertex (i);
				coCoord[i] = aaPos[vVrt[i]];
			}

		//	evaluate finite volume geometry
			geo.update (elem, coCoord, domain.subset_handler().get ());
			size_t noc = geo.num_scv ();

		//	loop corners to get the volumes and average the gradients
			for (size_t i = 0; i < noc; i++)
				aaVolume[vVrt[i]] += geo.scv(i).volume ();
		}
	}
	
#	ifdef UG_PARALLEL
	AttachmentAllReduce<Vertex> (grid, aVolume, PCL_RO_SUM);
#	endif
}

/**
 * computes the gradient at nodes of one grid element
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::compute_elem_grad
(
	DimFV1Geometry<dim> & geo, ///< FV geometry
	number uValue [], ///< values of the solution at the corners
	MathVector<dim> co_grad [], ///< for the gradients at the corners
	number * lsf, ///< corner values of the LSF (or NULL if none)
	CplUserData<number,dim> * if_val_data, ///< computes the values at the interface (if not NOLL)
	int si ///< subset index (only used if if_val_data != NULL)
)
{
	size_t noc = geo.num_scv ();
	
//	check the LSF
	if (lsf == NULL || lsf_sign (noc, lsf) != 0) // if no LSF or not intersected
	{
	//	loop corners to get the gradients
		for (size_t i = 0; i < noc; i++)
		{
		//	get scv
			const typename DimFV1Geometry<dim>::SCV& scv = geo.scv (i);

		//	sum up gradients of shape functions in corner
			co_grad[i] = 0.0;
			for (size_t sh = 0; sh < noc; sh++)
				VecScaleAppend (co_grad[i], uValue[sh], scv.global_grad (sh));
		}
	}
	else // the element is intersected, process every corner separately
	{
	//	loop corners to get the gradients
		for (size_t i = 0; i < noc; i++)
		{
		//	skip interface vertices
			if (lsf[i] <= lsf_threshold () && lsf[i] >= - lsf_threshold ())
			{
			//	this gradient is never really used: the value where it can occur are reset
				co_grad[i] = 0;
				continue;
			}
			
		//	extrapolate the solution
			number extValue [maxNumCo];
			extrapolate_by_lsf (if_val_data, si, geo, uValue, lsf, i, extValue);
			
		//	get scv
			const typename DimFV1Geometry<dim>::SCV& scv = geo.scv (i);

		//	sum up gradients of shape functions in corner
			co_grad[i] = 0.0;
			for (size_t sh = 0; sh < noc; sh++)
				VecScaleAppend (co_grad[i], extValue[sh], scv.global_grad (sh));
		}
	}
}

/**
 * computes the gradient at grid vertices and saves the result in the grid attachments
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::compute_vertex_grad
(
	TGridFunction& u, ///< gradient of this function
	DimFV1Geometry<dim>& geo, ///< FV geometry object
	t_aaVol& aaVolume, ///< the volumes
	ADimVector& aGradient, ///< where to save
	t_aaGrad& aaGradient, ///< accessor for the attachment
	CplUserData<number,dim> * if_val_data ///< computes the values at the interface (if not NOLL)
)
{
//	get domain
	domain_type& domain = *u.domain().get ();

//	get grid of domain
	grid_type& grid = *domain.grid ();

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor ();
	
//	local algebra
	std::vector<DoFIndex> multInd;
	LocalIndices locInd; LocalVector locU; LocalVector locLSF;

//	initialize attachment value
	SetAttachmentValues (aaGradient, grid.vertices_begin (), grid.vertices_end (), 0);

//	local values
	MathVector<dim> coCoord[maxNumCo];
	Vertex* vVrt[maxNumCo];
	number uValue[maxNumCo];
	MathVector<dim> globalGrad [maxNumCo];
	number lsfValue[maxNumCo], * lsf;

//	sum up all contributions of the sub control volumes to one vertex in an attachment
	for (int si = 0; si < domain.subset_handler()->num_subsets (); si++)
	{
	//TODO: Skipping boundary here can lead to wrong computation of the velocity at the boundary, cannot it?
	//	skip boundary
	//	if (m_dirichlet_sg.size () != 0) if (m_dirichlet_sg.contains (si)) continue;
	//	if (m_neumann_sg.size () != 0) if (m_neumann_sg.contains (si)) continue;
		
	//	loop grid elements of the full dimensionality
		ElemIterator iterEnd = u.template end<ElemType> (si);
		for (ElemIterator iter = u.template begin<ElemType> (si); iter != iterEnd; ++iter)
		{
		//	get Elem
			ElemType* elem = *iter;

		//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices ();
			for (size_t i = 0; i < numVertices; i++)
			{
				vVrt[i] = elem->vertex (i);
				coCoord[i] = aaPos[vVrt[i]];
			}

		//	evaluate finite volume geometry
			geo.update (elem, coCoord, domain.subset_handler().get ());
			size_t noc = geo.num_scv ();

		//	get the local solution
			u.indices (elem, locInd);
			locU.resize (locInd);
			GetLocalVector (locU, u);
			for (size_t i = 0; i < noc; i++)
				uValue[i] = locU (0, i);
			
		//	get the local LSF (if any)
			if (m_spLSF.valid ())
			{
				locLSF.resize (locInd);
				GetLocalVector (locLSF, *m_spLSF);
				for (size_t i = 0; i < noc; i++)
					lsfValue[i] = locLSF (0, i);
				lsf = lsfValue;
			}
			else
				lsf = NULL;

		//	get the gradient in the element
			compute_elem_grad (geo, uValue, globalGrad, lsf, if_val_data, si);

		//	loop corners to get the volumes and average the gradients
			for (size_t i = 0; i < noc; i++)
			{
			//	scale gradient by volume
				globalGrad[i] *= geo.scv(i).volume ();
			
			//	add it to the nodal gradient
				aaGradient[vVrt[i]] += globalGrad[i];
			}
		}
	}
	
#	ifdef UG_PARALLEL
	AttachmentAllReduce<Vertex> (grid, aGradient, PCL_RO_SUM);
#	endif

//	divide the gradients by the volumes
	for (int si = 0; si < u.num_subsets (); si++)
	{
	    for (VertexConstIterator iter = u.template begin<Vertex> (si);
								iter != u.template end<Vertex> (si); ++iter)
	    {
	    //	get vertex
		    Vertex* vrt = *iter;
		    if (aaVolume[vrt] != 0) //TODO: Eliminate this!
		        aaGradient[vrt] /= aaVolume[vrt];
	    }
	}
}

/**
 * sets Dirichlet values in solution vector for vertices in a given subset
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::assign_dirichlet
(
	TGridFunction& numsol
)
{
//	get domain of grid function
	domain_type& domain = *numsol.domain().get ();
	position_accessor_type aaPos = domain.position_accessor ();

//	loop the Dirichlet subsets
	std::vector<DoFIndex> ind (1);
	for (size_t i = 0; i < m_dirichlet_sg.size (); i++)
	{
		const int si = m_dirichlet_sg[i];
		
		for (VertexConstIterator iter = numsol.template begin<Vertex> (si);
									   iter != numsol.template end<Vertex> (si); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
			number exactVal;

		//	get vector holding all indices on the vertex
			numsol.inner_dof_indices (vrt, 0, ind);
		
		//	get the bc and save it in the solution
			if (m_imDirichlet.valid ())
				(*m_imDirichlet) (&exactVal, &aaPos[vrt], m_time, si, 1);
			else
				exactVal = DoFRef (*m_oldSol, ind[0]);
			DoFRef (numsol, ind[0]) = exactVal;
		}
	}
}

/**
 * slope limiter
 *
 * limit previously computed gradient so that the control-volume-wise linear
 * interpolation function does not introduce new maxima or minima into the data 
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::limit_grad
(
	TGridFunction& uOld, ///< old solution
	t_aaGrad& aaGrad ///< the gradient to limit
)
{
	grid_type& grid = *uOld.domain()->grid ();
	position_accessor_type& aaPos = uOld.domain()->position_accessor ();

	std::vector<DoFIndex> ind;

//	create Attachment for scv-volume size
	ANumber aMax;
	ANumber aMin;

//	attach to grid
	grid.attach_to_vertices (aMin);
	grid.attach_to_vertices (aMax);

//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaMin (grid, aMin);
	Grid::VertexAttachmentAccessor<ANumber> aaMax (grid, aMax);
	for (int si = 0; si < uOld.num_subsets(); si++)
	{
		for (VertexConstIterator iter = uOld.template begin<Vertex> (si);
				iter != uOld.template end<Vertex> (si); ++iter)
		{
			Vertex* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
			//	read indices on vertex
			//	get vector holding all indices on the vertex
			uOld.inner_dof_indices (vrt, 0, ind);
			aaMax[vrt] = DoFRef (uOld, ind[0]);
			aaMin[vrt] = DoFRef (uOld, ind[0]);
		}
	}
	for (int si = 0; si < uOld.num_subsets (); si++)
	{
	//	Loop the edges
		for (EdgeConstIterator iter = uOld.template begin<Edge> (si);
				iter != uOld.template end<Edge> (si); ++iter)
		{
			Edge* edge = *iter;
			Vertex* vi = edge->vertex (0);
			Vertex* vj = edge->vertex (1);
			uOld.inner_dof_indices (vi, 0, ind);
			number ui = DoFRef (uOld, ind[0]);
			uOld.inner_dof_indices (vj, 0, ind);
			number uj = DoFRef (uOld, ind[0]);
			//UG_LOG ("edge " << aaPos[vi] << "-" << aaPos[vj] << " [" << ui << " " << uj << "]\n");
			if (uj<aaMin[vi])
				aaMin[vi] = uj;
			if (uj>aaMax[vi])
				aaMax[vi] = uj;
			if (ui<aaMin[vj])
				aaMin[vj] = ui;
			if (ui>aaMax[vj])
				aaMax[vj] = ui;
		}
	}
	for (int si = 0; si < uOld.num_subsets (); si++) //TODO Why originally only for Subset 0?
	{
	//	skip the Dirichlet boundary
	//TODO: Is this true here?
		if (m_dirichlet_sg.size () != 0) if (m_dirichlet_sg.contains (si)) continue;
		
	//	Loop the edges
		for (EdgeConstIterator iter = uOld.template begin<Edge> (si);
				iter != uOld.template end<Edge> (si); ++iter)
		{
			Edge* edge = *iter;
			Vertex* vi = edge->vertex (0);
			Vertex* vj = edge->vertex (1);
			MathVector<dim> coordi,coordj,coordij,distVec,gradi,gradj;
			gradi = aaGrad[vi];
			gradj = aaGrad[vj];
			coordi = aaPos[vi];
			coordj = aaPos[vj];
			uOld.inner_dof_indices (vi, 0, ind);
			number ui = DoFRef (uOld, ind[0]);
			uOld.inner_dof_indices (vj, 0, ind);
			number uj = DoFRef (uOld, ind[0]);
			VecScaleAdd (coordij,0.5,coordi,0.5,coordj);
			VecSubtract (distVec, coordij,coordi);
			number uij = ui + distVec*gradi;
			number alpha = 1;
			if (uij > ui)
			{
				if (uij > aaMax[vi]) alpha = (aaMax[vi] - ui) / (distVec * gradi);
				if (alpha < 1)
				{
					//UG_LOG ("edge " << coordi << " " << coordj << "\n");
					//UG_LOG (coordi << " u " << ui << " uij "  << uij << " max " << aaMax[vi] << " alpha " << alpha << "\n");
					aaGrad[vi] *= alpha;
				}
			}
			else
			{
				if (uij < aaMin[vi]) alpha = (aaMin[vi] - ui) / (distVec * gradi);
				if (alpha < 1)
				{
					//UG_LOG ("edge " << coordi << " " << coordj << "\n");
					//UG_LOG (coordi << " u " << ui << " uij "  << uij << " min " << aaMax[vi] << " alpha " << alpha << "\n");
					aaGrad[vi] *= alpha;
				}
			}
			VecSubtract (distVec, coordij, coordj);
			uij = uj + distVec*gradj;
			alpha = 1;
			if (uij > uj)
			{
				if (uij > aaMax[vj]) alpha = (aaMax[vj] - uj) / (distVec * gradj);
				if (alpha < 1)
				{
					//UG_LOG ("-- edge " << coordi << " " << coordj << "\n");
					//UG_LOG (coordj << " u " << uj << " uij "  << uij << " max " << aaMax[vj] << " alpha " << alpha << "\n");
					aaGrad[vj] *= alpha;
				}
			}
			else
			{
				if (uij < aaMin[vj]) alpha = (aaMin[vj] - uj) / (distVec * gradj);
				if (alpha < 1)
				{
					//UG_LOG ("-- edge " << coordi << " " << coordj << "\n");
					//UG_LOG (coordj << " u " << uj << " uij "  << uij << " min " << aaMax[vj] << " alpha " << alpha << "\n");
					aaGrad[vj] *= alpha;
				}
			}
         // UG_LOG (" coord vertex 0 " << aaPos[v0] << " coord vertex 1 " << aaPos[v1] << "\n");
		}

	//	Loop full-dimensional grid elements
		Vertex* vVrt [maxNumCo];
		number u [maxNumCo];
		MathVector<dim> coCoord [maxNumCo];
		MathVector<dim> grad [maxNumCo];
		ElemIterator iterEnd = uOld.template end<ElemType> (si);
		for (ElemIterator iter = uOld.template begin<ElemType> (si); iter != iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;

			//	get position accessor
			const position_accessor_type& aaPos = uOld.domain()->position_accessor ();

			//	compute center of mass
			MathVector<dim> center;
			center=0;
			size_t noc = elem->num_vertices ();
			for (size_t i = 0; i < noc; i++)
			{
				vVrt[i] = elem->vertex (i);
				coCoord[i] = aaPos[vVrt[i]];
				grad[i] = aaGrad[vVrt[i]] ;
				VecAppend (center,coCoord[i]);
				uOld.inner_dof_indices (vVrt[i], 0, ind);
				u[i] = DoFRef (uOld, ind[0]);
			}
			center /= noc;
			for (size_t i = 0; i < noc; i++)
			{
				number alpha=1;
				MathVector<dim> distVec;
				number uCenter;
				VecSubtract (distVec, center, coCoord[i]);
				uCenter = u[i] + distVec * grad[i];
				if (uCenter > u[i])
				{
					if (uCenter > aaMax[vVrt[i]]) alpha = (aaMax[vVrt[i]] - u[i]) / (distVec * grad[i]);
					if (alpha < 1)
					{
						// UG_LOG ("* " << coCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " max " << aaMax[vVrt[i]] << " alpha " << alpha << "\n");
						aaGrad[vVrt[i]] *= alpha;
					}
				}
				else
				{
					if (uCenter < aaMin[vVrt[i]]) alpha = (aaMin[vVrt[i]] - u[i]) / (distVec * grad[i]);
					if (alpha < 1)
					{
						// UG_LOG ("*#* " << coCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " min " << aaMin[vVrt[i]] << " alpha " << alpha << "\n");
						aaGrad[vVrt[i]] *= alpha;
					}
				}
			}
		}
	}
	//	detach from grid
	grid.detach_from_vertices (aMin);
	grid.detach_from_vertices (aMax);
}

/**
 * Computes the time steps of the discretization of the level-set equation.
 *
 * Main function for assembling and solving ls equation as described in 
 * Frolkovic/Mikula, HIGH-RESOLUTION FLUX-BASED LEVEL SET METHOD, SIAM
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::advect ()
{
//	should we compute anything?
	if (m_nrOfSteps <= 0) return;
	
//	get the grid functions
	if (m_oldSol.invalid () || m_newSol.invalid ())
		UG_THROW ("Grid functions for the solutions not specified.");
	TGridFunction& uOld = *m_oldSol;
	TGridFunction& uNew = *m_newSol;
	
//	get domain of grid function
	domain_type& domain = *uNew.domain().get ();
	
//	get the position accessor
    position_accessor_type aaPos = domain.position_accessor ();
	
//	get grid of domain
	grid_type& grid = *domain.grid ();

//	attachment for scv-volume size
	ANumber aScvVolume;

//	attachment for gradient
	ADimVector aGradient; // gradient of the solution
	ADimVector aVelGrad; // gradient used for the velocity (if any)
	ADimVector aSDFGrad; // gradient of the signed-distance function (if any)
	
//	attachments for the updates of the solution (due to the avective terms and the source)
	ANumber aUpdate;
	ANumber aSrc;
	
//	attachment for the update of the SDF (if needed)
	ANumber aSDFUpdate;
	
//	attachment to mark corners of intersected elements
	ABool aCoIE; // Corner of Intersected Element

//	attach to grid
	grid.attach_to_vertices (aScvVolume);
	grid.attach_to_vertices (aGradient);
	grid.attach_to_vertices (aUpdate);
	grid.attach_to_vertices (aSrc);
	grid.attach_to_vertices (aCoIE);

//	get attachment accessor to access values
	t_aaVol aaVolume (grid, aScvVolume);
	t_aaGrad aaGradient (grid, aGradient);
	t_aaUpd aaUpdate (grid, aUpdate);
	t_aaUpd aaSrc (grid, aSrc);
	t_aaCoIE aaCoIE (grid, aCoIE);
	t_aaGrad aaVelGrad;
	t_aaGrad aaSDFGrad;
	t_aaUpd aaSDFUpdate;
	
//	FV geometry
	DimFV1Geometry<dim> geo;

//	specify the neumann (outflow) bnd subsets for the geometry, so that the
//	geometry produces boundary faces (BF) for all sides of the
//	element, that is in one of the subsets
	for (size_t i = 0; i < m_neumann_sg.size (); i++)
		geo.add_boundary_subset (m_neumann_sg[i]);

//	calculate scv volume
	compute_volumes (uNew, geo, aScvVolume, aaVolume);
	
//	attachment for the gradient used for the velocity
	if (m_spVelPot.valid () && m_spVelPot != m_oldSol)
	{
	//	attach, access and compute
		grid.attach_to_vertices (aVelGrad);
		aaVelGrad.access (grid, aVelGrad);
		compute_vertex_grad (*m_spVelPot, geo, aaVolume, aVelGrad, aaVelGrad);
		if (m_limiter)
			limit_grad (*m_spVelPot, aaVelGrad);
	}
	else
	{
	//	use the gradient of the old solution for it
		aaVelGrad.access (grid, aGradient); // merely redirect the accessor
	}
	
//	attachment for the gradient of the SDF
	if (m_spLSF.valid ())
	{
		if (m_spSDF.invalid ())
			UG_THROW ("Computation with the LSF interface is only possible with the SDF. Specify it!");
		
		if (m_spSDF == m_oldSol)
		{ // merely redirect the accessors
			aaSDFGrad.access (grid, aGradient);
			aaSDFUpdate.access (grid, aUpdate);
		}
		else if (m_spSDF == m_spVelPot)
		{
		// merely redirect the accessor for the gradient
			aaSDFGrad.access (grid, aVelGrad);
		//	attach and access the update
			grid.attach_to_vertices (aSDFUpdate);
			aaSDFUpdate.access (grid, aSDFUpdate);
		}
		else
		{
		//	attach, access and compute the gradient
			grid.attach_to_vertices (aSDFGrad);
			aaSDFGrad.access (grid, aSDFGrad);
			compute_vertex_grad (*m_spSDF, geo, aaVolume, aSDFGrad, aaSDFGrad);
			if (m_limiter)
				limit_grad (*m_spSDF, aaSDFGrad);
		//	attach and access the update
			grid.attach_to_vertices (aSDFUpdate);
			aaSDFUpdate.access (grid, aSDFUpdate);
		}
	}
	
//	mark the corners of the intersected elements
	mark_CoIE (grid, aCoIE, aaCoIE);

    MathVector<dim> coord;
	std::vector<DoFIndex> ind;

//	local indices and values
	LocalIndices locInd; LocalVector locOldU, locLSF, locVelPot;
	
	if (this->debug_writer_valid())
	{
		this->write_debug (uNew, "HR-FB-LSM-sol_step_0");
	}
	
//	computation of the time steps
	size_t step = 1;
	VecAssign (uOld, uNew);
	m_CFL = 0; // the Courant number over all the steps
	while (true)
	{
	//	set the Courant number for OUTPUT to 0
		number current_CFL = 0;
		if (m_spCourant.valid ())
		{
			(* m_spCourant) = 0.0;
			#ifdef UG_PARALLEL
			m_spCourant->set_storage_type (PST_ADDITIVE);
			#endif
		}
		else
		{
			if (m_time_control)
				UG_THROW ("Time control is only possible if a grid function is specified for CFL (by save_CourantNumber_to)");
		}
	
	//	indicators of wrong signes of the gradient
		bool wrong_sgn_at_if_A = false, wrong_sgn_at_if_B = false;
		
	//	compute scv volume and the gradient
	    compute_vertex_grad (uOld, geo, aaVolume, aGradient, aaGradient, m_imInterfaceVal.get ());
	    if (m_limiter)
	    	limit_grad (uOld, aaGradient);
	    
	//	initialize attachment values
		SetAttachmentValues (aaUpdate, grid.vertices_begin (), grid.vertices_end (), 0);
		SetAttachmentValues (aaSrc, grid.vertices_begin (), grid.vertices_end (), 0);
		if (m_spLSF.valid () && m_spSDF != m_oldSol)
			SetAttachmentValues (aaSDFUpdate, grid.vertices_begin (), grid.vertices_end (), 0);

	//	loop over subsets to compute the new solution
	    for (int si = 0; si < uOld.num_subsets (); si++)
	    {
		//	loop elements compute the update of the solution
		    ElemIterator iterEnd = uNew.template end<ElemType> (si);
		    for (ElemIterator iter = uNew.template begin<ElemType> (si); iter != iterEnd; ++iter)
		    {
		    	ElemType* elem = *iter;
		    	uNew.indices (elem, locInd);
		    	locOldU.resize (locInd);
		    	
		    	GetLocalVector (locOldU, uOld);
		    	if (m_spLSF.invalid ())
			    	assemble_element (elem, geo, domain, locOldU, aaGradient, aaVelGrad, aaVolume, 0, aaUpdate, &aaSrc);
			    else
			    {
			    	locLSF.resize (locInd); locVelPot.resize (locInd);
		    		int sign = assemble_cut_element
		    			(elem, geo, domain, locOldU, locLSF, locVelPot,
		    				aaGradient, aaVelGrad, aaVolume, aaUpdate, &aaSrc, m_imInterfaceVal.get (), si);
		    		if (sign != 0)
			    		assemble_element (elem, geo, domain, locOldU,
			    			aaGradient, aaVelGrad, aaVolume, sign, aaUpdate, &aaSrc);
			    }
			}
			
		//	loop elements to assemble the update of the SDF (if needed)
			if (m_spLSF.valid () && m_spSDF != m_oldSol)
				for (ElemIterator iter = uNew.template begin<ElemType> (si); iter != iterEnd; ++iter)
				{
		    		ElemType* elem = *iter;
		    		
				//	check whether we need this element: we assemble only the elements with marked corners
					for (size_t co = 0; co < elem->num_vertices (); co++)
					if (aaCoIE [elem->vertex (co)])
					{
						m_spSDF->indices (elem, locInd);
						locOldU.resize (locInd); locLSF.resize (locInd); locVelPot.resize (locInd);
			
						GetLocalVector (locOldU, *m_spSDF);
						int sign = assemble_cut_element
							(elem, geo, domain, locOldU, locLSF, locVelPot,
								aaSDFGrad, aaVelGrad, aaVolume, aaSDFUpdate, NULL, NULL, si);
						if (sign != 0)
							assemble_element (elem, geo, domain, locOldU,
								aaSDFGrad, aaVelGrad, aaVolume, sign, aaSDFUpdate, NULL);
						
						break;
					}
				}
	    }
	    
#		ifdef UG_PARALLEL
		AttachmentAllReduce<Vertex> (grid, aUpdate, PCL_RO_SUM);
		AttachmentAllReduce<Vertex> (grid, aSrc, PCL_RO_SUM);
		if (m_spLSF.valid () && m_spSDF != m_oldSol)
			AttachmentAllReduce<Vertex> (grid, aSDFUpdate, PCL_RO_SUM);
		if (m_spCourant.valid ())
			m_spCourant->change_storage_type (PST_CONSISTENT);
#		endif
	    
	//	check the CFL
		if (m_spCourant.valid ())
		{
			current_CFL = 0;
			for (VertexConstIterator iter = m_spCourant->template begin<Vertex> ();
					iter != m_spCourant->template end<Vertex> (); ++iter)
			{
				m_spCourant->inner_dof_indices (*iter, 0, ind);
				number local_CFL = DoFRef (*m_spCourant, ind[0]);
				if (local_CFL > current_CFL)
					current_CFL = local_CFL;
			}
#			ifdef UG_PARALLEL
			{
				pcl::ProcessCommunicator procComm;
				current_CFL = procComm.allreduce (current_CFL, PCL_RO_MAX);
			}
#			endif
		}
		if (m_time_control && current_CFL > 1e-10 && (current_CFL < m_minCFL || current_CFL > m_maxCFL))
		{
			number old_dt = m_dt;
			m_dt = m_dt / current_CFL * (m_minCFL + m_maxCFL) / 2;
			if (m_bVerbose)
				UG_LOG ("CFL " << current_CFL << " ... resetting time step: " << old_dt << " -> " << m_dt << "\n");
			continue; // recompute the step
		}
	    
	//	take into account the source at the vertices and update the solution
		for (int si = 0; si < uOld.num_subsets (); si++)
		{
		//	skip the Dirichlet boundary
			if (m_dirichlet_sg.size () != 0) if (m_dirichlet_sg.contains (si)) continue;
			
			for (VertexConstIterator iter = uNew.template begin<Vertex> (si);
					    			iter != uNew.template end<Vertex> (si); ++iter)
			{
				Vertex* vrt = *iter;
				uNew.inner_dof_indices (vrt, 0, ind);
				if (!aaCoIE[vrt])
				{//	just a normal vertex, not at the interface
					number co_source;
					if (m_spLSF.valid ())
						co_source = (DoFRef (*m_spLSF, ind[0]) >= 0)? m_source_pos : m_source_neg;
					else
						co_source = m_source_pos;
					if (m_antiderivSrc) co_source *= aaSrc[vrt];
					DoFRef (uNew, ind[0]) += m_dt * (aaUpdate[vrt] + co_source);
				}
				else
				{//	at the interface; try to change the time step
					number lsf_val = DoFRef (*m_spLSF, ind[0]);
					number sdf_update = aaSDFUpdate[vrt];
					number dt_eff = - DoFRef (*m_spSDF, ind[0]) / sdf_update; // note that this is "-" normalized gradient!
					if (dt_eff > m_dt)
						dt_eff = m_dt;
					if (lsf_val > lsf_threshold ())
					{
						if (sdf_update >= 0)
						{ // this should not happen
							wrong_sgn_at_if_A = true;
							DoFRef (uNew, ind[0]) += m_dt * m_source_pos; // only the source term
						}
						else
						{
							number co_source = m_source_pos;
							if (m_antiderivSrc) co_source *= aaSrc[vrt];
							DoFRef (uNew, ind[0]) += dt_eff * (aaUpdate[vrt] + co_source);
						}
					}
					else if (lsf_val < - lsf_threshold ())
					{
						if (sdf_update <= 0)
						{ // this should not happen
							wrong_sgn_at_if_B = true;
							DoFRef (uNew, ind[0]) += m_dt * m_source_neg; // only the source term
						}
						else
						{
							number co_source = m_source_neg;
							if (m_antiderivSrc) co_source *= aaSrc[vrt];
							DoFRef (uNew, ind[0]) += dt_eff * (aaUpdate[vrt] + co_source);
						}
					}
					else // we consider the vertex as lying directly at the interface
					{
						if (m_imInterfaceVal.invalid ())
							DoFRef (uNew, ind[0]) = 0;
						else
							(* m_imInterfaceVal) (DoFRef (uNew, ind[0]), aaPos[vrt], m_time, si, vrt);
					}
				}
			}
	    }
	    
	//	warnings
		if (wrong_sgn_at_if_A)
			UG_LOG ("advect: Wrong (non-negative) sign of the gradient occured!\n");
		if (wrong_sgn_at_if_B)
			UG_LOG ("advect: Wrong (non-positive) sign of the gradient occured!\n");
	    
	//	the new solution computed
	    m_time += m_dt;
	    
	//	set the Dirichlet values
	    assign_dirichlet (uNew);
	    
	//	update the total CFL
		if (current_CFL > m_CFL) m_CFL = current_CFL;
	    
	//	the time step is done
		if (m_bVerbose)
		{
			UG_LOG ("step length: " << m_dt << " (step # " << step << ")\n");
			UG_LOG ("time: " << m_time << "\n");
			UG_LOG ("CFL in step: " << current_CFL << "\n");
			
			UG_LOG ("advect: TS # " << step << ", time " << m_time << " (dt = " << m_dt << ")");
			if (m_spCourant.valid ())
				UG_LOG (", CFL = " << current_CFL);
			UG_LOG ("\n");
		}
	
		if (this->debug_writer_valid())
		{
			char debug_file_name [128];
			sprintf (debug_file_name, "HR-FB-LSM-sol_step_%d", (int) step);
			this->write_debug (uNew, debug_file_name);
		}
        
    //	is that all?
        if (step >= m_nrOfSteps)
        	break;
        else
        {
        	step++;
	   		VecAssign (uOld, uNew);
	   	}
	}
	
	if (! m_bVerbose)
	{
		UG_LOG ("advect: " << step << " step(s) done, last dt: " << m_dt);
		if (m_spCourant.valid ())
			UG_LOG (", max. CFL = " << m_CFL);
		UG_LOG ("\n");
	}
	
//	detach from grid
    if (m_spVelPot.valid () && m_spVelPot != m_oldSol)
		grid.detach_from_vertices (aVelGrad);
	if (m_spLSF.valid () && m_spSDF != m_oldSol)
	{
		grid.detach_from_vertices (aSDFUpdate);
		if (m_spSDF != m_spVelPot)
			grid.detach_from_vertices (aSDFGrad);
	}
	grid.detach_from_vertices (aCoIE);
	grid.detach_from_vertices (aSrc);
	grid.detach_from_vertices (aUpdate);
	grid.detach_from_vertices (aGradient);
	grid.detach_from_vertices (aScvVolume);
}

/**
 * compute the normal velocity using a user-data object
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::compute_normal_vel
(
	SmartPtr<CplUserData<MathVector<dim>, dim> > spVelField, ///< user data for the velocity vector field
	SmartPtr<TGridFunction> spNormVel ///< to save the normal velocity
)
{
//	we need the solution
	if (! m_newSol.valid ())
		UG_THROW ("Specify the level-set function!");
	
//	get domain
	domain_type& domain = * (m_newSol->domain().get ());

//	get grid of domain
	grid_type& grid = *domain.grid ();

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor ();
	
//	FV geometry
	DimFV1Geometry<dim> geo;

//	specify the neumann (outflow) bnd subsets for the geometry, so that the
//	geometry produces boundary faces (BF) for all sides of the
//	element, that is in one of the subsets
	for (size_t i = 0; i < m_neumann_sg.size (); i++)
		geo.add_boundary_subset (m_neumann_sg[i]);

//	get the total scv volumes
	ANumber aScvVolume;
	grid.attach_to_vertices (aScvVolume);
	t_aaVol aaVolume (grid, aScvVolume);
	compute_volumes (*m_newSol, geo, aScvVolume, aaVolume);
	
//	attach, access and compute the gradient
	ADimVector aGradient;
	grid.attach_to_vertices (aGradient);
	t_aaGrad aaGradient (grid, aGradient);
	compute_vertex_grad (*m_newSol, geo, aaVolume, aGradient, aaGradient);
	if (m_limiter)
		limit_grad (*m_newSol, aaGradient);
		
//	local indices
	std::vector<DoFIndex> ind;
	
//	from now on, aScvVolume will keep the volumes unter the interface
	SetAttachmentValues (aaVolume, grid.vertices_begin (), grid.vertices_end (), 0);
	
//	reset the data
	(*spNormVel) = 0.0;

//	loop over subsets to sum up the normal velocities
	for (int si = 0; si < m_newSol->num_subsets (); si++)
	{
		ElemIterator iterEnd = m_newSol->template end<ElemType> (si);
		for (ElemIterator iter = m_newSol->template begin<ElemType> (si); iter != iterEnd; ++iter)
		{
			ElemType* elem = *iter;
			
		//	check if this is a "positive" element
			number uValue[maxNumCo];
			Vertex* vVrt[maxNumCo];
			for (size_t i = 0; i < elem->num_vertices (); ++i)
			{
				Vertex* vrt = elem->vertex (i);
				m_newSol->inner_dof_indices ((vVrt[i] = vrt), 0, ind);
				uValue[i] = DoFRef (*m_newSol, ind[0]);
			}
			if (lsf_sign (elem->num_vertices (), uValue) > 0)
				continue; // we do not consider this element
			
		//	get vertices and extract corner coordinates
			MathVector<dim> coCoord[maxNumCo];
			for (size_t i = 0; i < elem->num_vertices (); ++i)
				coCoord[i] = aaPos[vVrt[i]];

		//	update fv geometry
			geo.update (elem, coCoord, domain.subset_handler().get ());
			size_t noc = geo.num_scv ();
			
		//	get the nodal velocity
			MathVector<dim> co_vel[maxNumCo];
			(*spVelField) (co_vel, geo.scv_global_ips (), m_time, si,
				elem, geo.corners (), geo.scv_local_ips (), noc, NULL);
		
		//	compute the normal velocity
			for (size_t i = 0; i < noc; i++)
			{
				Vertex* vrt = vVrt[i];
				number vol = geo.scv(i).volume ();
    	    	number gnorm = VecLength (aaGradient [vrt]);
				spNormVel->inner_dof_indices (vrt, 0, ind);
				DoFRef (*spNormVel, ind[0]) += (co_vel [i] * aaGradient [vrt]) * vol / gnorm;
				aaVolume[vrt] += vol;
			}
		}
	}
	
#	ifdef UG_PARALLEL
	AttachmentAllReduce<Vertex> (grid, aScvVolume, PCL_RO_SUM);
	spNormVel->set_storage_type (PST_ADDITIVE);
	spNormVel->change_storage_type (PST_CONSISTENT);
#	endif
	
//	loop over subsets to divide the normal velocities by the volumes
	for (int si = 0; si < m_newSol->num_subsets (); si++)
	{
		for (VertexConstIterator iter = m_newSol->template begin<Vertex> (si);
								iter != m_newSol->template end<Vertex> (si); ++iter)
		{
			Vertex* vrt = *iter;
			number vol = aaVolume[vrt];
			if (vol == 0)
				continue; // we are not under the interface!
			spNormVel->inner_dof_indices (vrt, 0, ind);
			DoFRef (*spNormVel, ind[0]) /= vol;
		}
	}
	
//	detach from the grid
	grid.detach_from_vertices (aGradient);
	grid.detach_from_vertices (aScvVolume);
}

/**
 * Appends a vertical vector the the normal velocity (using a user-data object)
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::append_vertical_to_normal_vel
(
	SmartPtr<CplUserData<number, dim> > spZVelField, ///< user data for the vertical velocity field
	SmartPtr<TGridFunction> spNormVel ///< the normal velocity to append to
)
{
//	we need the solution
	if (! m_newSol.valid ())
		UG_THROW ("Specify the level-set function!");
	
//	get domain
	domain_type& domain = * (m_newSol->domain().get ());

//	get grid of domain
	grid_type& grid = *domain.grid ();

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor ();
	
//	FV geometry
	DimFV1Geometry<dim> geo;

//	specify the neumann (outflow) bnd subsets for the geometry, so that the
//	geometry produces boundary faces (BF) for all sides of the
//	element, that is in one of the subsets
	for (size_t i = 0; i < m_neumann_sg.size (); i++)
		geo.add_boundary_subset (m_neumann_sg[i]);

//	Attachment for the contribution of the normal velocity
	ANumber aNVelocity;
	grid.attach_to_vertices (aNVelocity);
	t_aaVol aaNVel (grid, aNVelocity);
	SetAttachmentValues (aaNVel, grid.vertices_begin (), grid.vertices_end (), 0);
	
//	get the total scv volumes
	ANumber aScvVolume;
	grid.attach_to_vertices (aScvVolume);
	t_aaVol aaVolume (grid, aScvVolume);
	compute_volumes (*m_newSol, geo, aScvVolume, aaVolume);
	
//	attach, access and compute the gradient
	ADimVector aGradient;
	grid.attach_to_vertices (aGradient);
	t_aaGrad aaGradient (grid, aGradient);
	compute_vertex_grad (*m_newSol, geo, aaVolume, aGradient, aaGradient);
	if (m_limiter)
		limit_grad (*m_newSol, aaGradient);
		
//	local indices
	std::vector<DoFIndex> ind;
	
//	from now on, aScvVolume will keep the volumes unter the interface
	SetAttachmentValues (aaVolume, grid.vertices_begin (), grid.vertices_end (), 0);
	
//	loop over subsets to sum up the normal velocities
	for (int si = 0; si < m_newSol->num_subsets (); si++)
	{
		ElemIterator iterEnd = m_newSol->template end<ElemType> (si);
		for (ElemIterator iter = m_newSol->template begin<ElemType> (si); iter != iterEnd; ++iter)
		{
			ElemType* elem = *iter;
			
		//	check if this is a "positive" element
			number uValue[maxNumCo];
			Vertex* vVrt[maxNumCo];
			for (size_t i = 0; i < elem->num_vertices (); ++i)
			{
				Vertex* vrt = elem->vertex (i);
				m_newSol->inner_dof_indices ((vVrt[i] = vrt), 0, ind);
				uValue[i] = DoFRef (*m_newSol, ind[0]);
			}
			if (lsf_sign (elem->num_vertices (), uValue) > 0)
				continue; // we do not consider this element
			
		//	get vertices and extract corner coordinates
			MathVector<dim> coCoord[maxNumCo];
			for (size_t i = 0; i < elem->num_vertices (); ++i)
				coCoord[i] = aaPos[vVrt[i]];

		//	update fv geometry
			geo.update (elem, coCoord, domain.subset_handler().get ());
			size_t noc = geo.num_scv ();
			
		//	get the nodal velocity
			number co_nVel[maxNumCo];
			(*spZVelField) (co_nVel, geo.scv_global_ips (), m_time, si,
				elem, geo.corners (), geo.scv_local_ips (), noc, NULL);
		
		//	compute the normal velocity
			for (size_t i = 0; i < noc; i++)
			{
				Vertex* vrt = vVrt[i];
				number vol = geo.scv(i).volume ();
    	    	number gnorm = VecLength (aaGradient [vrt]);
				aaNVel [vrt] += (co_nVel [i] * aaGradient [vrt] [dim-1]) * vol / gnorm;
				aaVolume[vrt] += vol;
			}
		}
	}
	
#	ifdef UG_PARALLEL
	AttachmentAllReduce<Vertex> (grid, aScvVolume, PCL_RO_SUM);
	AttachmentAllReduce<Vertex> (grid, aNVelocity, PCL_RO_SUM);
#	endif
	
//	loop over subsets to divide the normal velocities by the volumes
	for (int si = 0; si < m_newSol->num_subsets (); si++)
	{
		for (VertexConstIterator iter = m_newSol->template begin<Vertex> (si);
								iter != m_newSol->template end<Vertex> (si); ++iter)
		{
			Vertex* vrt = *iter;
			number vol = aaVolume[vrt];
			if (vol == 0)
				continue; // we are not under the interface!
			spNormVel->inner_dof_indices (vrt, 0, ind);
			DoFRef (*spNormVel, ind[0]) += aaNVel [vrt] / vol;
		}
	}
	
//	detach from the grid
	grid.detach_from_vertices (aGradient);
	grid.detach_from_vertices (aScvVolume);
	grid.detach_from_vertices (aNVelocity);
}

/**
 * set subsets of the Dirichlet boundary
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::set_dirichlet_boundary
(
	const char* subsets
)
{
	if (m_dirichlet_sg.subset_handler().invalid ())
	{
		if (m_newSol.invalid ())
			UG_THROW ("Grid function for the solution not set.");
		m_dirichlet_sg.set_subset_handler (m_newSol->subset_handler ());
	}
	
	m_dirichlet_sg.add (TokenizeString (subsets));
}

/**
 * set subsets of the Outflow boundary
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::set_outflow_boundary
(
	const char* subsets
)
{
	if (m_neumann_sg.subset_handler().invalid ())
	{
		if (m_newSol.invalid ())
			UG_THROW ("Grid function for the solution not set.");
		m_neumann_sg.set_subset_handler (m_newSol->subset_handler ());
	}
	
	m_neumann_sg.add (TokenizeString (subsets));
}

/**
 * checks whether two level-set functions specify the same interface
 */
template<typename TGridFunction>
void HiResFluxBasedLSM<TGridFunction>::compare_lsf_with
(
	SmartPtr<TGridFunction> spLSF2, ///< the second level-set function
	number eps ///< tolerance
)
{
	if (m_newSol.invalid ())
		UG_THROW ("compare_lsf_with: Solution is not specified!\n");
	
	bool check_failed = false;
	std::vector<DoFIndex> ind;
	
//	get domain
	domain_type& domain = * (m_newSol->domain().get ());

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor ();
	
// 1. Compare the values at the vertices
	UG_LOG ("--- Checking solution, phase 1: Nodal values\n");
	for (int si = 0; si < domain.subset_handler()->num_subsets (); si++)
	{
		UG_LOG (" -- subset "<< si << "\n");
		for (VertexConstIterator iter = m_newSol->template begin<Vertex> (si);
								iter != m_newSol->template end<Vertex> (si); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
		//	get the index in the vector
			m_newSol->inner_dof_indices (vrt, 0, ind);
		//	compare the values
			number lsf_1 = DoFRef (*m_newSol, ind[0]);
			number lsf_2 = DoFRef (*spLSF2, ind[0]);
			if (lsf_1 * lsf_2 <= 0 && (lsf_1 != 0 || lsf_2 != 0))
			{
				check_failed = true;
				UG_LOG ("  > vertex " << vrt->grid_data_index() << " @ (" << aaPos[vrt][0]);
				for (size_t i = 1; i < dim; i++)
					UG_LOG (", " << aaPos[vrt][i]);
				UG_LOG ("): sol = " << lsf_1 << ", lsf = " << lsf_2 << "\n");
			}
		}
	}
	if (check_failed)
	{
		UG_LOG ("--- failed.\n");
		return;
	}
	UG_LOG ("--- passed.\n");
	
//	2. Compare the interface in the elements
	UG_LOG ("--- Checking solution, phase 2: Interface in the elements\n");
	for (int si = 0; si < domain.subset_handler()->num_subsets (); si++)
	{
		UG_LOG (" -- subset "<< si << "\n");
	//	loop grid elements of the full dimensionality
		ElemIterator iterEnd = m_newSol->template end<ElemType> (si);
		for (ElemIterator iter = m_newSol->template begin<ElemType> (si); iter != iterEnd; ++iter)
		{
			ElemType* elem = *iter;
			const size_t numVertices = elem->num_vertices ();
			
		//	local values
			Vertex* vVrt[maxNumCo];
			number lsf_1 [maxNumCo], lsf_2 [maxNumCo];

		//	get the vertices, extract the values of the functions
			for (size_t i = 0; i < numVertices; i++)
			{
				vVrt[i] = elem->vertex (i);
				m_newSol->inner_dof_indices (vVrt[i], 0, ind);
				lsf_1[i] = DoFRef (*m_newSol, ind[0]);
				lsf_2[i] = DoFRef (*spLSF2, ind[0]);
			}
		
		//	check if this is an interface element
			if (lsf_sign (numVertices, lsf_1) != 0)
				continue;
			
		//	check the position of the interface
			bool bad_if_elem = false;
			number param_1 [maxNumCo][maxNumCo];
			number param_2 [maxNumCo][maxNumCo];
			memset (param_1, 0, maxNumCo * maxNumCo * sizeof (number));
			memset (param_2, 0, maxNumCo * maxNumCo * sizeof (number));
			for (size_t i = 0; i < numVertices; i++)
				for (size_t j = i + 1; j < numVertices; j++)
				{
					if (lsf_1[i] * lsf_1[j] > 0) continue;
					param_1[i][j] = lsf_1[i] / (lsf_1[i] - lsf_1[j]);
					param_2[i][j] = lsf_2[i] / (lsf_2[i] - lsf_2[j]);
					if (std::fabs (param_1[i][j] - param_2[i][j]) >= eps)
						bad_if_elem = true;
				}
			if (bad_if_elem)
			{
				check_failed = true;
				UG_LOG ("  > elem " << elem->grid_data_index() << " (" << numVertices << " corners):\n");
				for (size_t i = 0; i < numVertices; i++)
				{
					UG_LOG ("  : " << i << " @ (" << aaPos[vVrt[i]][0]);
					for (size_t j = 1; j < dim; j++)
						UG_LOG (", " << aaPos[vVrt[i]][j]);
					UG_LOG ("): sol = " << lsf_1[i] << ", lsf = " << lsf_2[i] << "\n");
				}
				for (size_t i = 0; i < numVertices; i++)
					for (size_t j = i + 1; j < numVertices; j++)
						if (std::fabs (param_1[i][j] - param_2[i][j]) >= eps)
							UG_LOG ("  ! " << i << " to " << j << ": sol = " << param_1[i][j] << ", lsf = " << param_2[i][j] << "\n");
			}
		}
	}
	if (check_failed)
	{
		UG_LOG ("--- failed.\n");
		return;
	}
	UG_LOG ("--- passed.\n");
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */
