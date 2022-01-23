/*
 * Copyright (c) 2021:  G-CSC, Goethe University Frankfurt
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
 * Linkers for extracting data near the level set
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET_FILTER_LINKER_H__
#define __H__UG__PLUGINS__LEVEL_SET_FILTER_LINKER_H__

#include "common/common.h"

#include "lib_disc/spatial_disc/dom_disc_embb.h"
#include "lib_disc/spatial_disc/user_data/linker/linker.h"

namespace ug {
namespace LevelSet {

/**
 * Class for the linker that extracts data on the cut elements.
 *
 * This linker evaluates the associated object only on the elements cut by the level set.
 * for the other elements, it returns 0.
 */
template <typename TDomain, typename TAlgebra, typename TData>
class LSFilterLinker
:	public StdDataLinker< LSFilterLinker<TDomain, TAlgebra, TData>, TData, TDomain::dim>
{
	typedef StdDataLinker< LSFilterLinker<TDomain, TAlgebra, TData>, TData, TDomain::dim> base_type;
	
public:

	//	world dimension of grid function
	static const int dim = TDomain::dim;

	//	domain type
	typedef TDomain domain_type;
	
	//	algebra type
	typedef TAlgebra algebra_type;
	
	//	extrapolation type
	typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;
	
	//	type of the data
	typedef TData data_type;
	
	//	constructor
	LSFilterLinker
	(
		SmartPtr<CplUserData<data_type, dim> > spData, ///< the original data to extract from
		SmartPtr<extrapol_type> spExtrapol ///< extrapolation by the LSF
	)
	:	m_spExtrapolation (spExtrapol)
	{
		if (spData.invalid ())
			UG_THROW ("LSFilterLinker: No valid data specified.");
		this->set_num_input(1);
		m_spData = spData;
		m_spDData = spData.template cast_dynamic<DependentUserData<data_type, dim> > ();
		this->set_input (0, spData, spData);
	}

	inline void evaluate (data_type& value,
						  const MathVector<dim>& globIP,
						  number time, int si) const
	{
		UG_THROW ("LSFilterLinker needs the grid element.");
	}

	template <int refDim>
	inline void evaluate
	(
		data_type vValue[],
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
				(rRefElem.num(0), elem, si, false, vCornerCoords, time) == 0)
			{ // element is cut
				(* m_spData) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				return;
			}
		}
		
		// element is not cut, set the value to 0
		for (size_t ip = 0; ip < nip; ++ip)
			vValue[ip] = 0;
	}
	
	template <int refDim>
	void eval_and_deriv
	(
		data_type vValue[],
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
		std::vector<std::vector<data_type> > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		//	check if the element is inside/cut/outside
		if (m_spExtrapolation.valid ())
		{
			const ReferenceObjectID roid = elem->reference_object_id ();
			const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim> (roid);
			
			if (((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, false, vCornerCoords, time) == 0)
			{ // element is cut
				const data_type* vV = m_spData->values (s);
				for (size_t i = 0; i < nip; i++)
				{
					vValue[i] = vV[i];
			
					if (! bDeriv) continue;
					for (size_t fct = 0; fct < m_spDData->num_fct(); fct++)
					{
						const data_type* vDData = m_spDData->deriv (s, i, fct);
						const size_t c_fct = this->input_common_fct (0, fct);
						for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
							vvvDeriv[i][c_fct][sh] = vDData [sh];
					}
				}
				return;
			}
		}
		
		//	element not cut, set the output to zero
		for (size_t ip = 0; ip < nip; ++ip)
			vValue[ip] = 0;
		if (bDeriv && ! this->zero_derivative ())
			this->set_zero (vvvDeriv, nip);
	}
	
private:

	//	extrapolation by the level-set function
	SmartPtr<extrapol_type> m_spExtrapolation;

	//	the data for the cut elements
	SmartPtr<CplUserData<data_type, dim> > m_spData;
	SmartPtr<DependentUserData<data_type, dim> > m_spDData;
};

/**
 * Class for the linker that switches between data on the cut and uncut elements.
 *
 * This linker evaluates the associated object only on the elements cut by the level set
 * and the other object for the other elements.
 */
template <typename TDomain, typename TAlgebra, typename TData>
class LSFilterLinker2
:	public StdDataLinker< LSFilterLinker2<TDomain, TAlgebra, TData>, TData, TDomain::dim>
{
	typedef StdDataLinker< LSFilterLinker<TDomain, TAlgebra, TData>, TData, TDomain::dim> base_type;
	
public:

	//	world dimension of grid function
	static const int dim = TDomain::dim;

	//	domain type
	typedef TDomain domain_type;
	
	//	algebra type
	typedef TAlgebra algebra_type;
	
	//	extrapolation type
	typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;
	
	//	type of the data
	typedef TData data_type;
	
	//	constructor
	LSFilterLinker2
	(
		SmartPtr<CplUserData<data_type, dim> > spData_ce, ///< the original data to extract from for cut elements
		SmartPtr<CplUserData<data_type, dim> > spData_we, ///< the original data to extract from for whole elements
		SmartPtr<extrapol_type> spExtrapol ///< extrapolation by the LSF
	)
	:	m_spExtrapolation (spExtrapol)
	{
		if (spData_ce.invalid () || spData_we.invalid ())
			UG_THROW ("LSFilterLinker2: No valid data specified.");
		this->set_num_input(2);
		
		m_spData_ce = spData_ce;
		m_spDData_ce = spData_ce.template cast_dynamic<DependentUserData<data_type, dim> > ();
		this->set_input (0, spData_ce, spData_ce);
		
		m_spData_we = spData_we;
		m_spDData_we = spData_we.template cast_dynamic<DependentUserData<data_type, dim> > ();
		this->set_input (1, spData_we, spData_we);
	}

	inline void evaluate (data_type& value,
						  const MathVector<dim>& globIP,
						  number time, int si) const
	{
		UG_THROW ("LSFilterLinker2 needs the grid element.");
	}

	template <int refDim>
	inline void evaluate
	(
		data_type vValue[],
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
				(rRefElem.num(0), elem, si, false, vCornerCoords, time) == 0)
			{ // element is cut
				(* m_spData_ce) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
				return;
			}
		}
		
		// element is not cut, set the value to 0
		(* m_spData_we) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
	}
	
	template <int refDim>
	void eval_and_deriv
	(
		data_type vValue[],
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
		std::vector<std::vector<data_type> > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
		//	check if the element is inside/cut/outside
		if (m_spExtrapolation.valid ())
		{
			const ReferenceObjectID roid = elem->reference_object_id ();
			const DimReferenceElement<dim>& rRefElem = ReferenceElementProvider::get<dim> (roid);
			
			if (((extrapol_type *) m_spExtrapolation.get())->check_elem_lsf
				(rRefElem.num(0), elem, si, false, vCornerCoords, time) == 0)
			{ // element is cut
				const data_type* vV = m_spData_ce->values (s);
				for (size_t i = 0; i < nip; i++)
				{
					vValue[i] = vV[i];
			
					if (! bDeriv) continue;
					for (size_t fct = 0; fct < m_spDData_ce->num_fct(); fct++)
					{
						const data_type* vDData = m_spDData_ce->deriv (s, i, fct);
						const size_t c_fct = this->input_common_fct (0, fct);
						for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
							vvvDeriv[i][c_fct][sh] = vDData [sh];
					}
				}
				return;
			}
		}
		
		//	element not cut, use the second data set
		const data_type* vV = m_spData_we->values (s);
		for (size_t i = 0; i < nip; i++)
		{
			vValue[i] = vV[i];
	
			if (! bDeriv) continue;
			for (size_t fct = 0; fct < m_spDData_we->num_fct(); fct++)
			{
				const data_type* vDData = m_spDData_we->deriv (s, i, fct);
				const size_t c_fct = this->input_common_fct (0, fct);
				for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
					vvvDeriv[i][c_fct][sh] = vDData [sh];
			}
		}
	}
	
private:

	//	extrapolation by the level-set function
	SmartPtr<extrapol_type> m_spExtrapolation;

	//	the data for the cut elements
	SmartPtr<CplUserData<data_type, dim> > m_spData_ce;
	SmartPtr<DependentUserData<data_type, dim> > m_spDData_ce;

	//	the data for the whole elements
	SmartPtr<CplUserData<data_type, dim> > m_spData_we;
	SmartPtr<DependentUserData<data_type, dim> > m_spDData_we;
};

} // namespace LevelSet
} // end namespace ug

#endif // __H__UG__PLUGINS__LEVEL_SET_FILTER_LINKER_H__

/* End of File */
