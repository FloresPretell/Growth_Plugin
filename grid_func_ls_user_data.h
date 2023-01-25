/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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
 * grid_func_ls_user_data.h
 */

#ifndef __H__UG__PLUGINS__GRID_FUNC_LS_USER_DATA__
#define __H__UG__PLUGINS__GRID_FUNC_LS_USER_DATA__

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/spatial_disc/dom_disc_embb.h"

namespace ug{
namespace LevelSet {

/**
 * GridFunction data under the level set.
 *
 * This is the extension of the GridFunctionNumberData for the Level-Set method:
 * This class sets the values over the level set to 0.
 */
template <typename TGridFunction>
class GridFuncLSNumberData
: public StdDependentUserData<GridFuncLSNumberData<TGridFunction>,
  number, TGridFunction::dim>
{
	public:
		//	world dimension of grid function
		static const int dim = TGridFunction::dim;
		
		//	domain type
		typedef typename TGridFunction::domain_type domain_type;
		
		//	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;
		
		//	extrapolation type
		typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

		private:
		//	grid function
		SmartPtr<TGridFunction> m_spGridFct;
		
		//	extrapolation by the level-set function
		SmartPtr<extrapol_type> m_spExtrapolation;

		//	component of function
		size_t m_fct;

		//	local finite element id
		LFEID m_lfeID;

	public:
		/// constructor
		GridFuncLSNumberData(SmartPtr<TGridFunction> spGridFct, const char* cmp, SmartPtr<extrapol_type> spExtrapol)
		: m_spGridFct(spGridFct), m_spExtrapolation(spExtrapol)
		{
			this->set_functions(cmp);

			//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

			//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("GridFuncLSNumberData: Function space does not contain"
						" a function with name " << cmp << ".");

			//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);
			
			if (m_lfeID != LFEID (LFEID::LAGRANGE, dim, 1))
				UG_THROW ("GridFuncLSNumberData: Only Lagrange spaces are currently supported.");
		};

		virtual bool continuous() const
		{
			return LocalFiniteElementProvider::continuous(m_lfeID);
		}

		template <int refDim>
		void eval_and_deriv(number vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<number> > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();
			const DimReferenceElement<refDim>& rRefElem = ReferenceElementProvider::get<refDim>(roid);

			//	check if the element is inside/cut/outside
			int inside = 1;
			if(m_spExtrapolation.valid())
			{
				int g_level = m_spGridFct->grid_level().level();
				if(g_level == GridLevel::TOP)
					g_level = m_spGridFct->approx_space()->num_levels() - 1;
				if((inside = ((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
					(rRefElem.num(0), elem, si, g_level, false, vCornerCoords, time)) < 0)
				{
					for(size_t ip = 0; ip < nip; ++ip) // element outside
						vValue[ip] = 0; //TODO: set the interface value
					return;
				}
			}
			
			//	get trial space
			try{
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
						LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

				//	memory for shapes
				std::vector<number> vShape;
				std::vector<number> vVal;

				//	get multiindices of element
				std::vector<DoFIndex> ind;
				m_spGridFct->dof_indices(elem, m_fct, ind);

				//	loop ips
				for(size_t ip = 0; ip < nip; ++ip)
				{
					//	evaluate at shapes at ip
					rTrialSpace.shapes(vShape, vLocIP[ip]);
					
					//	get the values at the corners and extrapolate them
					vVal.resize(vShape.size());
					for(size_t sh = 0; sh < vVal.size(); ++sh)
						vVal[sh] = DoFRef(*m_spGridFct, ind[sh]);
					if(inside == 0) // the element is cut
						m_spExtrapolation->extrapolate_by_lsf (rRefElem.num(0), &(vVal[0]), m_fct);

					// 	compute solution at integration point
					vValue[ip] = 0.0;
					for(size_t sh = 0; sh < vShape.size(); ++sh)
						vValue[ip] += vVal[sh] * vShape[sh];
				}

				if(bDeriv){
					for(size_t ip = 0; ip < nip; ++ip){
						//	evaluate at shapes at ip
						rTrialSpace.shapes(vShape, vLocIP[ip]);

						for(size_t sh = 0; sh < vShape.size(); ++sh)
							vvvDeriv[ip][0][sh] = vShape[sh];
					}
				}
			}
			UG_CATCH_THROW("GridFuncLSNumberData: Shape Function Set missing for"
					" Reference Object: "<<roid<<", Trial Space: "
					<<m_lfeID<<", refDim="<<refDim);

		}
};

/**
 * GridFunction gradient data under the level set.
 *
 * This is the extension of the GridFunctionGradientData for the Level-Set method:
 * This class sets the gradient over the level set to 0.
 */
template <typename TGridFunction>
class GridFuncLSGradientData
: public StdDependentUserData<GridFuncLSGradientData<TGridFunction> ,
  MathVector<TGridFunction::dim>, TGridFunction::dim>
{
	public:
	//	world dimension of grid function
	static const int dim = TGridFunction::dim;

	//	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	//	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;
	
	//	extrapolation type
	typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

	private:
	// grid function
	SmartPtr<TGridFunction> m_spGridFct;

	//	extrapolation by the level-set function
	SmartPtr<extrapol_type> m_spExtrapolation;

	//	component of function
	size_t m_fct;

	//	local finite element id
	LFEID m_lfeID;

	public:
	/// constructor
	GridFuncLSGradientData(SmartPtr<TGridFunction> spGridFct, const char* cmp, SmartPtr<extrapol_type> spExtrapol)
	: m_spGridFct(spGridFct), m_spExtrapolation(spExtrapol)
	{
		this->set_functions(cmp);

		//	get function id of name
		m_fct = spGridFct->fct_id_by_name(cmp);

		//	check that function exists
		if(m_fct >= spGridFct->num_fct())
			UG_THROW("GridFuncLSGradientData: Function space does not contain"
					" a function with name " << cmp << ".");

		//	local finite element id
		m_lfeID = spGridFct->local_finite_element_id(m_fct);
		
		if (m_lfeID != LFEID (LFEID::LAGRANGE, dim, 1))
			UG_THROW ("GridFuncLSGradientData: Only Lagrange spaces are currently supported.");
	};

	virtual bool continuous() const
	{
		return false;
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
		const ReferenceObjectID roid = elem->reference_object_id();
		const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim>(roid);
		
		//	check if the element is inside/cut/outside
		int inside = 1;
		if(m_spExtrapolation.valid())
		{
			int g_level = m_spGridFct->grid_level().level();
			if(g_level == GridLevel::TOP)
				g_level = m_spGridFct->approx_space()->num_levels() - 1;
			if((inside = ((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, g_level, false, vCornerCoords, time)) < 0)
			{
				for(size_t ip = 0; ip < nip; ++ip) // element outside
					VecSet(vValue[ip], 0);
				return;
			}
		}
		
		//	get trial space
		try{
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

			//	storage for shape function at ip
			std::vector<MathVector<refDim> > vLocGrad;
			MathVector<refDim> locGrad;
			
			//	storage for corner values
			std::vector<number> vVal;

			//	Reference Mapping
			MathMatrix<dim, refDim> JTInv;
			
			//	get multiindices of element
			std::vector<DoFIndex > ind;
			m_spGridFct->dof_indices(elem, m_fct, ind);
	
			// treat the cut elements separately
			if(inside != 0) // the element is NOT cut
			{
				//	get reference element mapping by reference object id
				std::vector<MathMatrix<refDim, dim> > vJTTmp(nip);
				if(vJT == NULL)
				{
					DimReferenceMapping<refDim, dim>& mapping
						= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);

					//	compute transformation matrices
					mapping.jacobian_transposed(&(vJTTmp[0]), vLocIP, nip);

					//	store tmp Gradient
					vJT = &(vJTTmp[0]);
				}

				//	loop ips
				for(size_t ip = 0; ip < nip; ++ip)
				{
					//	evaluate at shapes at ip
					rTrialSpace.grads(vLocGrad, vLocIP[ip]);
	
					//	get the values at the corners
					vVal.resize(vLocGrad.size());
					for(size_t sh = 0; sh < vVal.size(); ++sh)
						vVal[sh] = DoFRef(*m_spGridFct, ind[sh]);
	
					//	compute grad at ip
					VecSet(locGrad, 0.0);
					for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
						VecScaleAppend(locGrad, vVal[sh], vLocGrad[sh]);
	
					Inverse(JTInv, vJT[ip]);
					MatVecMult(vValue[ip], JTInv, locGrad);
				}
			}
			else // REMARK: All the ip's get the same gradient!
			{	
				const DimReferenceElement<refDim>& rRefElem
					= ReferenceElementProvider::get<refDim> (roid);
				MathVector<dim> baseGrad, elemGrad;
				MathMatrix<refDim, dim> JT;
				
				DimReferenceMapping<refDim, dim>& mapping
						= ReferenceMappingProvider::get<refDim, dim>(roid, vCornerCoords);
				
				VecSet(elemGrad, 0.0);
				
				//	loop base corners
				size_t n_base_co = 0;
				for(size_t base_co = 0; base_co < rRefElem.num(0); ++base_co)
				if(m_spExtrapolation->corner_inside (base_co))
				{
					//	evaluate at shapes at ip
					rTrialSpace.grads(vLocGrad, rRefElem.corner(base_co));
	
					//	get the values at the corners and extrapolate them
					vVal.resize(vLocGrad.size());
					for(size_t sh = 0; sh < vVal.size(); ++sh)
						vVal[sh] = DoFRef(*m_spGridFct, ind[sh]);
					m_spExtrapolation->extrapolate_by_lsf (rRefElem.num(0), base_co, &(vVal[0]), m_fct);
	
					//	compute the local gradient at ip
					VecSet(locGrad, 0.0);
					for(size_t sh = 0; sh < vLocGrad.size(); ++sh)
						VecScaleAppend(locGrad, vVal[sh], vLocGrad[sh]);
					
					//	compute the gradient
					mapping.jacobian_transposed(JT, rRefElem.corner(base_co));
					Inverse(JTInv, JT);
					MatVecMult(baseGrad, JTInv, locGrad);
					
					//	add the contribution to the element gradient
					elemGrad += baseGrad;
					
					++n_base_co;
				}
				if(n_base_co == 0)
					UG_THROW("GridFunctionNumberData: Failed in a cut element.");
				
				//	average the gradient and copy it to all the ip's
				elemGrad /= (number) n_base_co;
				for(size_t ip = 0; ip < nip; ++ip)
					vValue[ip] = elemGrad;

			}
		}
		UG_CATCH_THROW("GridFunctionNumberData: Failed for"
				" Reference Object: "<<roid<<", Trial Space: "
				<<m_lfeID<<", refDim="<<refDim);

		if(bDeriv)
			UG_THROW("Not implemented.");
	}
};

/**
 * Extension of the Bear-Scheidegger dispersion for the Level-Set method.
 *
 * This class sets the dispersion over AND AT the interface to 0.
 */
template <typename TDomain, typename TAlgebra>
class LSBearScheidegger
	: public StdDataLinker< LSBearScheidegger<TDomain, TAlgebra>, MathMatrix<TDomain::dim,TDomain::dim>, TDomain::dim>
{
	typedef StdDataLinker< LSBearScheidegger<TDomain, TAlgebra>, MathMatrix<TDomain::dim,TDomain::dim>, TDomain::dim> base_type;
	
public:

	//	world dimension of grid function
	static const int dim = TDomain::dim;

	//	domain type
	typedef TDomain domain_type;
	
	//	algebra type
	typedef TAlgebra algebra_type;
	
	//	extrapolation type
	typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;
	
	//	constructor
	LSBearScheidegger
	(
		SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > spDispersion, ///< dispersion
		SmartPtr<extrapol_type> spExtrapol ///< extrapolation by the LSF
	)
	:	m_spExtrapolation (spExtrapol)
	{
		if (spDispersion.invalid ())
			UG_THROW ("LSBearScheidegger: No valid disperion specified.");
		this->set_num_input(1);
		m_spDispersion = spDispersion;
		m_spDDispersion = spDispersion.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> > ();
		this->set_input (0, spDispersion, spDispersion);
	}

	inline void evaluate (MathMatrix<dim,dim>& D,
						  const MathVector<dim>& globIP,
						  number time, int si) const
	{
		UG_THROW ("LSBearScheidegger needs the grid element.");
	}

	template <int refDim>
	inline void evaluate
	(
		MathMatrix<dim,dim> vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		//	check if the element is inside/cut/outside
		if (m_spExtrapolation.valid ())
		{
			const ReferenceObjectID roid = elem->reference_object_id ();
			const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim> (roid);
			
			if (((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, false, vCornerCoords, time) <= 0)
			{
				for (size_t ip = 0; ip < nip; ++ip) // element is cut or outside
					MatSet(vValue[ip], 0);
				return;
			}
		}
		
		//	compute the dispersion by the usual formula
		(* m_spDispersion) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
	}
	
	template <int refDim>
	void eval_and_deriv
	(
		MathMatrix<dim,dim> vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		bool bDeriv,
		int s,
		std::vector<std::vector<MathMatrix<dim,dim> > > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		//	check if the element is inside/cut/outside
		if (m_spExtrapolation.valid ())
		{
			const ReferenceObjectID roid = elem->reference_object_id ();
			const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim> (roid);
			
			if (((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, false, vCornerCoords, time) <= 0)
			{
				for (size_t ip = 0; ip < nip; ++ip) // element is cut or outside
					MatSet(vValue[ip], 0);
				if (bDeriv && ! this->zero_derivative ())
					this->set_zero (vvvDeriv, nip);
				return;
			}
		}
		
		//	compute the dispersion and its derivatives by the usual formulas
		const MathMatrix<dim,dim>* vVector = m_spDispersion->values (s);
		for (size_t i = 0; i < nip; i++)
		{
			vValue[i] = vVector[i];
			
			if (! bDeriv) continue;
			for (size_t fct = 0; fct < m_spDDispersion->num_fct(); fct++)
			{
				const MathMatrix<dim,dim>* vDDispersion = m_spDDispersion->deriv (s, i, fct);
				const size_t c_fct = this->input_common_fct (0, fct);
				for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
					vvvDeriv[i][c_fct][sh] = vDDispersion [sh];
			}
		}
	}
	
private:

	//	extrapolation by the level-set function
	SmartPtr<extrapol_type> m_spExtrapolation;

	//	the original dispersion
	SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > m_spDispersion;
	SmartPtr<DependentUserData<MathMatrix<dim,dim>, dim> > m_spDDispersion;
};

/**
 * Extension of the Bear-Scheidegger dispersion for the Level-Set method.
 *
 * This class sets the dispersion over the interface to 0. At the interface, it sets
 * the different dispersion as under the interface
 */
template <typename TDomain, typename TAlgebra>
class LSBearScheidegger2
	: public StdDataLinker< LSBearScheidegger2<TDomain, TAlgebra>, MathMatrix<TDomain::dim,TDomain::dim>, TDomain::dim>
{
	typedef StdDataLinker< LSBearScheidegger2<TDomain, TAlgebra>, MathMatrix<TDomain::dim,TDomain::dim>, TDomain::dim> base_type;
	
public:

	//	world dimension of grid function
	static const int dim = TDomain::dim;

	//	domain type
	typedef TDomain domain_type;
	
	//	algebra type
	typedef TAlgebra algebra_type;
	
	//	extrapolation type
	typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;
	
	//	constructor
	LSBearScheidegger2
	(
		SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > spDispersion_under, ///< dispersion under the interface
		SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > spDispersion_at, ///< dispersion at the interface
		SmartPtr<extrapol_type> spExtrapol ///< extrapolation by the LSF
	)
	:	m_spExtrapolation (spExtrapol)
	{
		if (spDispersion_under.invalid ())
			UG_THROW ("LSBearScheidegger2: No valid disperion under the interface specified.");
		if (spDispersion_at.invalid ())
			UG_THROW ("LSBearScheidegger2: No valid disperion at the interface specified.");
		this->set_num_input(2);
		m_spDispersion_under = spDispersion_under;
		m_spDDispersion_under = spDispersion_under.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> > ();
		this->set_input (0, spDispersion_under, spDispersion_under);
		m_spDispersion_at = spDispersion_at;
		m_spDDispersion_at = spDispersion_at.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> > ();
		this->set_input (1, spDispersion_at, spDispersion_at);
	}

	inline void evaluate (MathMatrix<dim,dim>& D,
						  const MathVector<dim>& globIP,
						  number time, int si) const
	{
		UG_THROW ("LSBearScheidegger2 needs the grid element.");
	}

	template <int refDim>
	inline void evaluate
	(
		MathMatrix<dim,dim> vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		//	check if the element is inside/cut/outside
		if (m_spExtrapolation.valid ())
		{
			const ReferenceObjectID roid = elem->reference_object_id ();
			const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim> (roid);
			int lsf;
			
			if ((lsf = ((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, false, vCornerCoords, time)) < 0)
			{
				for (size_t ip = 0; ip < nip; ++ip) // element is cut outside
					MatSet(vValue[ip], 0);
				return;
			}
			else if (lsf == 0)
			{
				(* m_spDispersion_at) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				return;
			}
		}
		
		//	compute the dispersion by the usual formula
		(* m_spDispersion_under) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
	}
	
	template <int refDim>
	void eval_and_deriv
	(
		MathMatrix<dim,dim> vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		bool bDeriv,
		int s,
		std::vector<std::vector<MathMatrix<dim,dim> > > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		//	check if the element is inside/cut/outside
		if (m_spExtrapolation.valid ())
		{
			const ReferenceObjectID roid = elem->reference_object_id ();
			const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim> (roid);
			int lsf;
			
			if ((lsf = ((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, false, vCornerCoords, time)) < 0)
			{
				for (size_t ip = 0; ip < nip; ++ip) // element is cut outside
					MatSet(vValue[ip], 0);
				if (bDeriv && ! this->zero_derivative ())
					this->set_zero (vvvDeriv, nip);
				return;
			}
			else if (lsf == 0)
			{
				//	compute the dispersion and its derivatives by the usual formulas at the interface
				const MathMatrix<dim,dim>* vVector = m_spDispersion_at->values (s);
				for (size_t i = 0; i < nip; i++)
				{
					vValue[i] = vVector[i];
			
					if (! bDeriv) continue;
					for (size_t fct = 0; fct < m_spDDispersion_at->num_fct(); fct++)
					{
						const MathMatrix<dim,dim>* vDDispersion_at = m_spDDispersion_at->deriv (s, i, fct);
						const size_t c_fct = this->input_common_fct (0, fct);
						for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
							vvvDeriv[i][c_fct][sh] = vDDispersion_at [sh];
					}
				}
				return;
			}
		}
		
		//	compute the dispersion and its derivatives by the usual formulas under the interface
		const MathMatrix<dim,dim>* vVector = m_spDispersion_under->values (s);
		for (size_t i = 0; i < nip; i++)
		{
			vValue[i] = vVector[i];
			
			if (! bDeriv) continue;
			for (size_t fct = 0; fct < m_spDDispersion_under->num_fct(); fct++)
			{
				const MathMatrix<dim,dim>* vDDispersion_under = m_spDDispersion_under->deriv (s, i, fct);
				const size_t c_fct = this->input_common_fct (0, fct);
				for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
					vvvDeriv[i][c_fct][sh] = vDDispersion_under [sh];
			}
		}
	}
	
private:

	//	extrapolation by the level-set function
	SmartPtr<extrapol_type> m_spExtrapolation;

	//	the original dispersion (under the interface)
	SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > m_spDispersion_under;
	SmartPtr<DependentUserData<MathMatrix<dim,dim>, dim> > m_spDDispersion_under;

	//	the effective dispersion at the interface
	SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > m_spDispersion_at;
	SmartPtr<DependentUserData<MathMatrix<dim,dim>, dim> > m_spDDispersion_at;
};

} // namespace LevelSet
} // end namespace ug

#endif /* __H__UG__PLUGINS__GRID_FUNC_LS_USER_DATA__ */
