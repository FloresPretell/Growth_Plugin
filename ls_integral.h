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
#ifndef __H__UG__PLUGINS__LEVEL_SET_LS_INTEGRAL_H__
#define __H__UG__PLUGINS__LEVEL_SET_LS_INTEGRAL_H__

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
 * This class computes the integral of a given grid function over the given full-dim. subsets.
 *
 * \tparam TGridFunc	type of the grid function with the values to integrade
 */
template <typename TGridFunction>
class LSIntegral
{
public:
///	type of the grid function
	typedef TGridFunction gf_type;
	
///	world dimension
	static const int dim = gf_type::dim;
	
///	this type
	typedef LSIntegral<gf_type> this_type;

///	domain type
	typedef typename gf_type::domain_type domain_type;
	
/// type of the position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	
///	algebra type
	typedef typename gf_type::algebra_type algebra_type;
	
///	algebra type for the LSF
	typedef CPUAlgebra lsf_algebra_type;
	
///	grid function type for the LSF
	typedef GridFunction<domain_type, lsf_algebra_type> ls_gf_type;
	
public:

///	Constructor
	LSIntegral
	(
		SmartPtr<gf_type>& sp_gf, ///< the grid function
		const char * fct_name, ///< 'function' to take the data from
		SmartPtr<ls_gf_type>& sp_lsf ///< the leve-set function (if any)
	)
	:	m_sp_gf (sp_gf), m_spLSF (sp_lsf)
	{
		if ((m_fct = sp_gf->fct_id_by_name (fct_name)) >= sp_gf->num_fct ())
			UG_THROW ("LSAveData: Function space does not contain any function with name '" << fct_name << "'.");
	
		if (sp_gf->local_finite_element_id (m_fct) != LFEID(LFEID::LAGRANGE, dim, 1))
			UG_THROW ("LSAveData: Only vertex-centered grid functions are supported.");
	}

///	Destructor
	virtual ~LSIntegral () {}
	
///	sets the subsets to restrict the computation on
	void on_subsets
	(
		const char * ss_names ///< names of the subsets
	)
	{
		m_ssGrp.set_subset_handler (m_sp_gf->domain()->subset_handler ());
		m_ssGrp.add (TokenizeString (ss_names));
	}

///	computes the integrals
	void compute ();
	
///	returns the integral over the negative part of subsets
	number integral_over_subsets
	(
		const char * ss_names ///< name of the subset
	) const;
	
///	return the integral over the entire negative subdomain
	number integral () const
	{
		return m_integral;
	}
	
private:

///	adds contributions of all elements of a given type
	template <typename TElem>
	void add_integrals_of_all ();
	
///	helper class for the computation of the volumes
	struct AddIntegrals
	{
		this_type * m_pThis;
		
	///	constructor
		AddIntegrals
		(
			this_type * pThis ///< pointer to the master class
		)
		: m_pThis (pThis) {}
		
	///	computation of the volume for all elements of the given type
		template <typename TElem> void operator() (TElem)
		{
			m_pThis->template add_integrals_of_all<TElem> ();
		}
	};
	
private:

	SmartPtr<gf_type> m_sp_gf; ///< the grid function
	size_t m_fct; ///< the function inside of the grid function
	
	SmartPtr<ls_gf_type> m_spLSF; ///< the specified LSF
	SubsetGroup m_ssGrp; ///< the subset group (if subsets specified)
	
	std::vector<number> m_ss_integral; ///< integrals over the subsets
	number m_integral; ///< the entire integral
};

} // namespace LevelSet
} // end namespace ug

#include "ls_integral_impl.h"

#endif // __H__UG__PLUGINS__LEVEL_SET_LS_INTEGRAL_H__

/* End of File */
