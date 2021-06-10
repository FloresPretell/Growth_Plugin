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

/**
 * Computation of an integral over a domain with the embedded boundary: implementation.
 */

#include "lib_disc/reference_element/reference_element_traits.h"
#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif

#include "ls_volume.h"

namespace ug{
namespace LevelSet{

/*---- Class 'LSIntegral': ----*/

/**
 * Sums up the (partial) integrals of all the elements of all the types.
 */
template <typename TGridFunc>
void LSIntegral<TGridFunc>::compute_for
(
	SmartPtr<gf_type> sp_gf, ///< the grid function
	const char * fct_name ///< 'function' to take the data from
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<dim>::DimElemList ElemList;
	
	m_sp_gf = sp_gf;
	if ((m_fct = sp_gf->fct_id_by_name (fct_name)) >= sp_gf->num_fct ())
		UG_THROW ("LSIntegral: Function space does not contain any function with name '" << fct_name << "'.");

	if (sp_gf->local_finite_element_id (m_fct) != LFEID(LFEID::LAGRANGE, dim, 1))
		UG_THROW ("LSIntegral: Only vertex-centered grid functions are supported.");
	
	int n_ss = m_spLSF->num_subsets ();

//	sum up all the integrals
	m_integral = 0;
	m_ss_integral.resize (n_ss);
	for (int si = 0; si < n_ss; si++) m_ss_integral[si] = 0;
	boost::mpl::for_each<ElemList> (AddIntegrals (this));
	
#ifdef UG_PARALLEL
//	sum up the volumes from different processes
	pcl::ProcessCommunicator procComm;
	m_integral = procComm.allreduce (m_integral, PCL_RO_SUM);
	
	for (int si = 0; si < n_ss; si++)
		m_ss_integral[si] = procComm.allreduce (m_ss_integral[si], PCL_RO_SUM);
#endif
}

/**
 * Extracts the integral over the given subsets
 */
template <typename TGridFunc>
number LSIntegral<TGridFunc>::integral_over_subsets
(
	const char * ss_names ///< name of the subset
) const
{
	SubsetGroup ss_grp (m_spLSF->domain()->subset_handler ());
	ss_grp.add (TokenizeString (ss_names));
	
	number ss_integral = 0;
	for (size_t i = 0; i < ss_grp.size (); i++)
		ss_integral += m_ss_integral[ss_grp[i]];
	return ss_integral;
}

/**
 * Sums up the (partial) integrals of all the elements of one type.
 */
template <typename TGridFunc>
template <typename TElem>
void LSIntegral<TGridFunc>::add_integrals_of_all ()
{
	typedef typename gf_type::template traits<TElem>::const_iterator ElemIter;
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_t;

	static const size_t num_corners = ref_elem_t::numCorners;
	
	const ls_gf_type & lsf = * m_spLSF;
	const gf_type & integrand = * m_sp_gf;
	const position_accessor_type & aaPos = lsf.domain()->position_accessor ();
	std::vector<DoFIndex> ind (1);
	MathVector<dim> corners [num_corners];
	number lsf_values [num_corners];
	
	for (int si = 0; si < lsf.num_subsets (); si++)
	{
		if (m_ssGrp.subset_handler().valid () && ! m_ssGrp.contains (si))
			continue; // skip this subset: it is not mentioned in the specified list
		
		number ss_integral = 0;
		ElemIter iterEnd = lsf.template end<TElem> (si);
		for (ElemIter iter = lsf.template begin<TElem> (si); iter != iterEnd; ++iter)
		{
			TElem * elem = *iter;
			
		//	get the corner coordinates ans the values of the LSF
			for (size_t i = 0; i < num_corners; i++)
			{
				Vertex * vrt = elem->vertex (i);
				corners [i] = aaPos [vrt];
				if (lsf.inner_dof_indices (vrt, 0, ind) != 1)
					UG_THROW ("LSIntegral: Not a scalar grid function for the LSF!");
				lsf_values [i] = DoFRef (lsf, ind [0]);
			}
			
		//	compute the volumes
			number vol_plus, vol_minus;
			if (LSElementSize<ref_elem_t, dim>::compute (corners, lsf_values, vol_plus, vol_minus) > 0)
				continue;
		
		//	get the average of the integrand
			number ave_integrand = 0;
			size_t n_co_minus = 0;
			for (size_t i = 0; i < num_corners; i++)
			if (lsf_values [i] < 0)
			{
				Vertex * vrt = elem->vertex (i);
				if (integrand.inner_dof_indices (vrt, m_fct, ind) != 1)
					UG_THROW ("LSIntegral: Not a scalar grid function for the integrand!");
				ave_integrand += DoFRef (integrand, ind[0]);
				n_co_minus++;
			}
			if (n_co_minus == 0) // this can happen because of the treatment of pos/neg nodes in LSElementSize
				continue;
			
		//	add the contribution to the subset
			ss_integral += vol_minus * ave_integrand / n_co_minus;
		}
		m_ss_integral[si] += ss_integral;
		m_integral += ss_integral;
	}
}

/*---- Class 'LSHeavisideIntegral': ----*/

/**
 * Sums up the (partial) integrals of all the elements of all the types.
 */
template <typename TGridFunc>
void LSHeavisideIntegral<TGridFunc>::compute_for
(
	SmartPtr<gf_type> sp_gf, ///< the grid function
	const char * fct_name, ///< 'function' to take the data from
	number limval ///< the value for the Heaviside function
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<dim>::DimElemList ElemList;
	
	m_sp_gf = sp_gf; m_limval = limval;
	if ((m_fct = sp_gf->fct_id_by_name (fct_name)) >= sp_gf->num_fct ())
		UG_THROW ("LSHeavisideIntegral: Function space does not contain any function with name '" << fct_name << "'.");

	if (sp_gf->local_finite_element_id (m_fct) != LFEID(LFEID::LAGRANGE, dim, 1))
		UG_THROW ("LSHeavisideIntegral: Only vertex-centered grid functions are supported.");
	
	int n_ss = m_spLSF->num_subsets ();

//	sum up all the integrals
	m_levol = m_gevol = 0;
	m_ss_levol.resize (n_ss); m_ss_gevol.resize (n_ss);
	for (int si = 0; si < n_ss; si++) m_ss_levol[si] = m_ss_gevol[si] = 0;
	boost::mpl::for_each<ElemList> (AddIntegrals (this));
	
#ifdef UG_PARALLEL
//	sum up the volumes from different processes
	pcl::ProcessCommunicator procComm;
	m_levol = procComm.allreduce (m_levol, PCL_RO_SUM);
	m_gevol = procComm.allreduce (m_gevol, PCL_RO_SUM);
	
	for (int si = 0; si < n_ss; si++)
	{
		m_ss_levol[si] = procComm.allreduce (m_ss_levol[si], PCL_RO_SUM);
		m_ss_gevol[si] = procComm.allreduce (m_ss_gevol[si], PCL_RO_SUM);
	}
#endif
}

/**
 * Extracts the volume for gf <= limval over the given subsets
 */
template <typename TGridFunc>
number LSHeavisideIntegral<TGridFunc>::le_vol_over_subsets
(
	const char * ss_names ///< name of the subset
) const
{
	SubsetGroup ss_grp (m_spLSF->domain()->subset_handler ());
	ss_grp.add (TokenizeString (ss_names));
	
	number ss_levol = 0;
	for (size_t i = 0; i < ss_grp.size (); i++)
		ss_levol += m_ss_levol[ss_grp[i]];
	return ss_levol;
}

/**
 * Extracts the volume for gf >= limval over the given subsets
 */
template <typename TGridFunc>
number LSHeavisideIntegral<TGridFunc>::ge_vol_over_subsets
(
	const char * ss_names ///< name of the subset
) const
{
	SubsetGroup ss_grp (m_spLSF->domain()->subset_handler ());
	ss_grp.add (TokenizeString (ss_names));
	
	number ss_gevol = 0;
	for (size_t i = 0; i < ss_grp.size (); i++)
		ss_gevol += m_ss_gevol[ss_grp[i]];
	return ss_gevol;
}

/**
 * Sums up the (partial) integrals of all the elements of one type.
 */
template <typename TGridFunc>
template <typename TElem>
void LSHeavisideIntegral<TGridFunc>::add_integrals_of_all ()
{
	typedef typename gf_type::template traits<TElem>::const_iterator ElemIter;
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_t;

	static const size_t num_corners = ref_elem_t::numCorners;
	
	const ls_gf_type & lsf = * m_spLSF;
	const gf_type & integrand = * m_sp_gf;
	const position_accessor_type & aaPos = lsf.domain()->position_accessor ();
	std::vector<DoFIndex> ind (1);
	MathVector<dim> corners [num_corners];
	number lsf_values [num_corners];
	
	for (int si = 0; si < lsf.num_subsets (); si++)
	{
		if (m_ssGrp.subset_handler().valid () && ! m_ssGrp.contains (si))
			continue; // skip this subset: it is not mentioned in the specified list
		
		number ss_levol = 0, ss_gevol = 0;
		ElemIter iterEnd = lsf.template end<TElem> (si);
		for (ElemIter iter = lsf.template begin<TElem> (si); iter != iterEnd; ++iter)
		{
			TElem * elem = *iter;
			
		//	get the corner coordinates ans the values of the LSF
			for (size_t i = 0; i < num_corners; i++)
			{
				Vertex * vrt = elem->vertex (i);
				corners [i] = aaPos [vrt];
				if (lsf.inner_dof_indices (vrt, 0, ind) != 1)
					UG_THROW ("LSHeavisideIntegral: Not a scalar grid function for the LSF!");
				lsf_values [i] = DoFRef (lsf, ind [0]);
			}
			
		//	compute the volumes
			number vol_plus, vol_minus;
			if (LSElementSize<ref_elem_t, dim>::compute (corners, lsf_values, vol_plus, vol_minus) > 0)
				continue;
		
		//	get the averaged integrand
			size_t sum_integrand = 0;
			size_t n_co_minus = 0;
			for (size_t i = 0; i < num_corners; i++)
			if (lsf_values [i] < 0)
			{
				Vertex * vrt = elem->vertex (i);
				if (integrand.inner_dof_indices (vrt, m_fct, ind) != 1)
					UG_THROW ("LSHeavisideIntegral: Not a scalar grid function for the integrand!");
				number gf_val = DoFRef (integrand, ind[0]);
				if (gf_val <= m_limval)
					sum_integrand += 1;
				n_co_minus++;
			}
			if (n_co_minus == 0) // this can happen because of the treatment of pos/neg nodes in LSElementSize
				continue;
			
		//	add the contribution to the subset
			ss_levol += (vol_minus * sum_integrand) / n_co_minus;
			ss_gevol += (vol_minus * (n_co_minus - sum_integrand)) / n_co_minus;
		}
		m_ss_levol[si] += ss_levol; m_ss_gevol[si] += ss_gevol;
		m_levol += ss_levol; m_gevol += ss_gevol;
	}
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */
