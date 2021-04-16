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
 * Tools for averaging of variable values in the elements (mainly for plotting)
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET_AVE_DATA_H__
#define __H__UG__PLUGINS__LEVEL_SET_AVE_DATA_H__

#include <vector>
#include <iostream>
#include <limits>
#include <cmath>

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/dom_disc_embb.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/element_list_traits.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"

namespace ug {
namespace LevelSet {

/**
 * Class for averaging of a grid function in the elements
 *
 * This class computes the average of a given grid function over the element.
 *
 * \tparam TGridFunc	type of the grid function with the values to averate
 */
template <typename TGridFunction>
class LSAveData
:	public StdDependentUserData<LSAveData<TGridFunction>, number, TGridFunction::dim>
{
public:
///	type of the grid function
	typedef TGridFunction gf_type;
	
///	world dimension
	static const int dim = gf_type::dim;
	
///	this type
	typedef LSAveData<gf_type> this_type;

///	the base type
	typedef StdDependentUserData<this_type, number, dim> base_type;
	
///	domain type
	typedef typename gf_type::domain_type domain_type;
	
///	algebra type
	typedef typename gf_type::algebra_type algebra_type;
	
///	extrapolation type
	typedef IInterfaceExtrapolation<domain_type, algebra_type> extrapol_type;

public:

///	Constructor
	LSAveData
	(
		SmartPtr<gf_type> sp_gf, ///< the grid function
		const char * fct_name ///< 'function' to take the data from
	)
	:	m_sp_gf (sp_gf)
	{
		if ((m_fct = sp_gf->fct_id_by_name (fct_name)) >= sp_gf->num_fct ())
			UG_THROW ("LSAveData: Function space does not contain any function with name '" << fct_name << "'.");
	
		if (sp_gf->local_finite_element_id (m_fct) != LFEID(LFEID::LAGRANGE, dim, 1))
			UG_THROW ("LSAveData: Only vertex-centered grid functions are supported.");
	}

///	Destructor
	virtual ~LSAveData () {}
	
///	Set the extrapolation (if any)
	void set_extrapolation (SmartPtr<extrapol_type> sp_extrapol) {m_sp_extrapol = sp_extrapol;}
		
///	These are element-based data, they are not continuous
	virtual bool continuous () const {return false;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return true;}

/// Performs the main computations:
	template <int refDim>
	void eval_and_deriv
	(
		number vValue[], ///< for the computed value
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
		std::vector<std::vector<number> > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	)
	{
	//	reference object id
		const ReferenceObjectID roid = elem->reference_object_id();
		const DimReferenceElement<dim>& r_ref_elem = ReferenceElementProvider::get<dim>(roid);
	
	//	check if the element is inside/cut/outside
		int inside = 1;
		if(m_sp_extrapol.valid ())
		{
			int g_level = m_sp_gf->grid_level().level();
			
			if (g_level == GridLevel::TOP)
				g_level = m_sp_gf->approx_space()->num_levels () - 1;
			
			if ((inside = ((extrapol_type *) m_sp_extrapol.get())->check_elem_lsf
				(r_ref_elem.num (0), elem, si, g_level, false, vCornerCoords, time)) < 0)
			{
				for(size_t ip = 0; ip < nip; ++ip) // element outside
					vValue[ip] = 0; //TODO: set the interface value
				return;
			}
		}
	
	//	get multiindices of element
		std::vector<DoFIndex> ind;
		m_sp_gf->dof_indices (elem, m_fct, ind);
		if (ind.size () != r_ref_elem.num (0))
			UG_THROW ("LSAveData: Wrong number of values per element.");

	//	average the values
		size_t n_vals = 0;
		number sum = 0;
		for (size_t co = 0; co < ind.size (); co++)
		{
			if (inside == 0 && ! m_sp_extrapol->corner_inside (co))
				continue;
			sum += DoFRef (*m_sp_gf, ind[co]);
			n_vals++;
		}
		if (n_vals != 0) // if the element is inside
			sum /= n_vals;
		else
			sum = 0; //TODO: set the interface value (this point should be unreachable)
	
	//	return the result
		for (size_t ip = 0; ip < nip; ip++)
			vValue[ip] = sum;
	}
	
private:

	SmartPtr<gf_type> m_sp_gf; ///< the grid function
	size_t m_fct; ///< the function inside of the grid function
	SmartPtr<extrapol_type> m_sp_extrapol; ///< the extrapolation (if any)
};

} // namespace LevelSet
} // end namespace ug

#endif // __H__UG__PLUGINS__LEVEL_SET_POSITION_H__

/* End of File */
