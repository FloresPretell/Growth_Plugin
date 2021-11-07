/*
 * Author: Dmitry Logashenko
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
 * Utilities for the computation of the vector velocity from the normal velocity.
 */

#ifndef __H__UG__PLUGINS__LEVEL_SET__NORMVEL_UTIL_H__
#define __H__UG__PLUGINS__LEVEL_SET__NORMVEL_UTIL_H__

// ug4 headers
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"

namespace ug{
namespace LevelSet{

/**
 * UserData based class that computes the velocity vector as the normalized gradient of a given function.
 */
template <typename TGridFunc>
class EikonalVel
	: public StdDependentUserData
				<EikonalVel<TGridFunc>, MathVector<TGridFunc::dim>, TGridFunc::dim>
{
public:
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

private:
/// 'potential' of the velocity: direction is its normalized gradient; typically the LSF
	SmartPtr<TGridFunc> m_spVelPot;

///	local finite element id (assumed to be the same for both the functions)
	LFEID m_lfeID;

public:
/// constructor
	EikonalVel
	(
		SmartPtr<TGridFunc> spPot ///< grid function with the 'potential' of the velocity
	)
	: m_spVelPot (spPot)
	{
	//	local finite element ids
		m_lfeID = m_spVelPot->local_finite_element_id(0);
	};
	
///	The vector field of the gradients is not continuous over the sides of the elements
	virtual bool continuous () const {return false;}
	
///	No derivatives implemented
	virtual bool zero_derivative () const {return true;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return true;}

/// Performs the main computations:
	template <int refDim>
	void eval_and_deriv
	(
		MathVector<dim> vValue[], ///< for the computed value
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject * elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector * u,
		bool bDeriv,
		int s,
		std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		const MathMatrix<refDim, dim> * vJT = NULL
	) const
	{
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("EikonalVel: Derivatives are not implemented.");
		
	//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();

	//	compute the local -> global trannsformation if it is not given
		std::vector<MathMatrix<refDim, dim> > vJTTmp(nip);
		if(vJT == NULL)
		{
			DimReferenceMapping<refDim, dim>& mapping
				= ReferenceMappingProvider::get<refDim, dim> (roid, vCornerCoords);
			mapping.jacobian_transposed(&(vJTTmp[0]), vLocIP, nip);
			vJT = &(vJTTmp[0]);
		}

	//	the shape functions
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

	//	memory for gradients and indices
		std::vector<DoFIndex> ind;
		std::vector<MathVector<refDim> > vLocGrad;
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		
	//	loop the integration points
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	compute the direction of the velocity
			MathVector<dim> grad;
			
		//	a) evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	b) get multiindices of element
			std::vector<DoFIndex > ind;
			m_spVelPot->dof_indices(elem, 0, ind);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				VecScaleAppend(locGrad, DoFRef (*m_spVelPot, ind[sh]), vLocGrad[sh]);

			RightInverse (JTInv, vJT[ip]);
			MatVecMult(grad, JTInv, locGrad);
		
		//	compute the velocity vector
			number len = VecLength (grad);
			if (len > 1e-15)
				VecScale (vValue[ip], grad, 1 / len);
			else
				VecSet (vValue[ip], 0.0);
		}
	};
};

/**
 * UserData based class that computes the vector velocity from the normal velocity.
 */
template <typename TGridFunc>
class VelByNormalVel
	: public StdDependentUserData
				<VelByNormalVel<TGridFunc>, MathVector<TGridFunc::dim>, TGridFunc::dim>
{
public:
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

private:
///	normal velocity
	SmartPtr<TGridFunc> m_spNV;

/// 'potential' of the velocity: direction is its normalized gradient; typically the LSF
	SmartPtr<TGridFunc> m_spVelPot;

///	local finite element id (assumed to be the same for both the functions)
	LFEID m_lfeID;

public:
/// constructor
	VelByNormalVel
	(
		SmartPtr<TGridFunc> spNV, ///< grid function with the normal velocity
		SmartPtr<TGridFunc> spPot ///< grid function with the 'potential' of the velocity
	)
	: m_spNV (spNV), m_spVelPot (spPot)
	{
	//	local finite element ids
		LFEID lfeID_nv = m_spNV->local_finite_element_id(0);
		LFEID lfeID_pot = m_spVelPot->local_finite_element_id(0);
		
		if (lfeID_nv != lfeID_pot)
			UG_THROW ("VelByNormalVel: Both grid functions should have the same approx. spaces");
		m_lfeID = lfeID_nv;
	};
	
///	The vector field of the gradients is not continuous over the sides of the elements
	virtual bool continuous () const {return false;}

///	No derivatives implemented
	virtual bool zero_derivative () const {return true;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return true;}

/// Performs the main computations:
	template <int refDim>
	void eval_and_deriv
	(
		MathVector<dim> vValue[], ///< for the computed value
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject * elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector * u,
		bool bDeriv,
		int s,
		std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		const MathMatrix<refDim, dim> * vJT = NULL
	) const
	{
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("VelByNormalVel: Derivatives are not implemented.");
		
	//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();

	//	compute the local -> global trannsformation if it is not given
		std::vector<MathMatrix<refDim, dim> > vJTTmp(nip);
		if(vJT == NULL)
		{
			DimReferenceMapping<refDim, dim>& mapping
				= ReferenceMappingProvider::get<refDim, dim> (roid, vCornerCoords);
			mapping.jacobian_transposed(&(vJTTmp[0]), vLocIP, nip);
			vJT = &(vJTTmp[0]);
		}

	//	the shape functions
		const LocalShapeFunctionSet<refDim>& rTrialSpace =
				LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

	//	memory for shapes, gradients and indices
		std::vector<number> vShape;
		std::vector<DoFIndex> ind;
		std::vector<MathVector<refDim> > vLocGrad;
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		
	//	loop the integration points
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	compute the normal velocity:
			number nv = 0;
		
		//	a) evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	b) get multiindices of element
			m_spNV->dof_indices(elem, 0, ind);

		// 	c) compute solution at integration point
			for(size_t sh = 0; sh < vShape.size(); ++sh)
				nv += DoFRef(*m_spNV, ind[sh]) * vShape[sh];
		
		//	compute the direction of the velocity
			MathVector<dim> grad;
			
		//	a) evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	b) get multiindices of element
			m_spVelPot->dof_indices(elem, 0, ind);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
				VecScaleAppend(locGrad, DoFRef (*m_spVelPot, ind[sh]), vLocGrad[sh]);

			RightInverse (JTInv, vJT[ip]);
			MatVecMult(grad, JTInv, locGrad);
		
		//	compute the velocity vector
			number len = VecLength (grad);
			if (len > 1e-15)
				VecScale (vValue[ip], grad, nv / len);
			else
				VecSet (vValue[ip], 0.0);
		}
	};
};

} // end namespace LevelSet
} // end namespace ug

#endif // __H__UG__PLUGINS__LEVEL_SET__NORMVEL_UTIL_H__

/* End of File */
