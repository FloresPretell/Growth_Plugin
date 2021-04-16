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

/**
 * Measuring volumes of the parts separaterd by the level set.
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET__LS_VOLUME_H__
#define __H__UG__PLUGINS__LEVEL_SET__LS_VOLUME_H__

// ug4 headers
#include "common/common.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#include "lib_disc/common/geometry_util.h"

namespace ug{
namespace LevelSet{

/**
 * A class for computation of the volume of the parts of the domain
 * separated by the level set.
 *
 * \tparam TGridFunc	type of the grid function for the LSF
 */
template <typename TGridFunc>
class LSVolume
{
public:

///	'this' type
	typedef LSVolume<TGridFunc> this_type;

///	grid function type
	typedef TGridFunc grid_func_type;

///	domain type
	typedef typename grid_func_type::domain_type domain_type;
	
///	algebra type
	typedef typename grid_func_type::algebra_type algebra_type;

/// type of the position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	
///	world dimension
	static const int dim = domain_type::dim;
	
public:

///	class constructor
	LSVolume
	(
		SmartPtr<grid_func_type> spLSF ///< the grid function with the LSF
	)
	:	m_spLSF (spLSF),
		m_volume_plus (-1), m_volume_minus (-1), // dummy values
		m_check_positivity (false)
	{};
	
///	sets the subsets to restrict the computation on
	void on_subsets
	(
		const char * ss_names ///< names of the subsets
	);

///	set the check-positivity flag
	void check_positivity (bool v) {m_check_positivity = v;}
	
///	computes the volumes
	void compute ();
	
///	extracts the volumes enclosed in given subsets
	void volume_in_subsets
	(
		const char * ss_names, ///< name of the subset
		number & vol_plus, ///< the "positive" volume
		number & vol_minus ///< the "negative" volume
	) const;
	
///	returns the volume in the positive part
	number volume_plus () const {return m_volume_plus;}
	
///	returns the volume in the positive part of subsets
	number volume_plus_in_subsets
	(
		const char * ss_names ///< name of the subset
	) const
	{
		number vol_plus, vol_minus;
		volume_in_subsets (ss_names, vol_plus, vol_minus);
		return vol_plus;
	}
	
/// returns the volume of the negative part
	number volume_minus () const {return m_volume_minus;}
	
///	returns the volume in the negative part of subsets
	number volume_minus_in_subsets
	(
		const char * ss_names ///< name of the subset
	) const
	{
		number vol_plus, vol_minus;
		volume_in_subsets (ss_names, vol_plus, vol_minus);
		return vol_minus;
	}
	
///	prints the details
	void print_details () const;
	
private:

///	adds contributions of all elements of a given type
	template <typename TElem>
	void add_volumes_of_all ();
	
///	helper class for the computation of the volumes
	struct AddVolumes
	{
		this_type * m_pThis;
		
	///	constructor
		AddVolumes
		(
			this_type * pThis ///< pointer to the master class
		)
		: m_pThis (pThis) {}
		
	///	computation of the volume for all elements of the given type
		template <typename TElem> void operator() (TElem)
		{
			m_pThis->template add_volumes_of_all<TElem> ();
		}
	};
	
private:

	SmartPtr<TGridFunc> m_spLSF; ///< the specified LSF
	
	SubsetGroup m_ssGrp; ///< the subset group (if subsets specified)
	
	number m_volume_plus; ///< computed volume in the "positive part" of the domain
	number m_volume_minus; ///< computed volume in the "negative part" of the domain
	
	std::vector<number> m_ss_vol_plus; ///< volume in the 'positive' subsets
	std::vector<number> m_ss_vol_minus; ///< volume in the 'negative' subsets
	
	bool m_check_positivity; ///< if to check the positivity of the volumes in the elements
};

/**
 * A class for the computation of the volumes of the parts of an element
 * intersected by the level-set. Note that the level set is considered as
 * planar in the element.
 *
 * \tparam TRefElem		type of the reference element
 * \tparam WDim			dimension of the world
 */
template <typename TRefElem, int WDim>
class LSElementSize
{
public:
	
	typedef TRefElem ref_element_type; ///< type of the reference element
	static const int dim = WDim; ///< dimension of the world
	
public:

///	computation of the volumes
	/**
	 * Return value: positive or negative if the entire element is in the corresponding
	 * subdomain; zero otherwise.
	 */
	static int compute
	(
		const MathVector<WDim> * corner, ///< coordinates of the corners
		const number * lsf, ///< values of the LSF at the corners
		number & vol_plus, ///< volume of the 'positive' part of the element
		number & vol_minus ///< volume of the 'negative' part of the element
	)
	{
		UG_THROW ("LSElementSize: Generic version not implemented");
		return 0; // never reached
	}
};

/**
 * Specialization of LSElementSize for edges
 * \see LSElementSize
 */
template <int WDim>
class LSElementSize<ReferenceEdge, WDim>
{
public:
	
	typedef ReferenceEdge ref_element_type; ///< type of the reference element
	static const int dim = WDim; ///< dimension of the world
	
public:

///	computation of the volumes (for an edge)
	/**
	 * Return value: positive or negative if the entire element is in the corresponding
	 * subdomain; zero otherwise.
	 */
	static int compute
	(
		const MathVector<WDim> * corner, ///< coordinates of the corners
		const number * lsf, ///< values of the LSF at the corners
		number & vol_plus, ///< volume of the 'positive' part of the element
		number & vol_minus ///< volume of the 'negative' part of the element
	);
};

/**
 * Specialization of LSElementSize for triangles
 * \see LSElementSize
 */
template <int WDim>
class LSElementSize<ReferenceTriangle, WDim>
{
public:
	
	typedef ReferenceTriangle ref_element_type; ///< type of the reference element
	static const int dim = WDim; ///< dimension of the world
	
public:

///	computation of the volumes (for a triangle)
	/**
	 * Return value: positive or negative if the entire element is in the corresponding
	 * subdomain; zero otherwise.
	 */
	static int compute
	(
		const MathVector<WDim> * corner, ///< coordinates of the corners
		const number * lsf, ///< values of the LSF at the corners
		number & vol_plus, ///< volume of the 'positive' part of the element
		number & vol_minus ///< volume of the 'negative' part of the element
	);
};

/**
 * Specialization of LSElementSize for tetrahedra
 * \see LSElementSize
 */
template <int WDim>
class LSElementSize<ReferenceTetrahedron, WDim>
{
public:
	
	typedef ReferenceTetrahedron ref_element_type; ///< type of the reference element
	static const int dim = WDim; ///< dimension of the world
	
public:

///	computation of the volumes (for a tetrahedrod)
	/**
	 * Return value: positive or negative if the entire element is in the corresponding
	 * subdomain; zero otherwise.
	 */
	static int compute
	(
		const MathVector<WDim> * corner, ///< coordinates of the corners
		const number * lsf, ///< values of the LSF at the corners
		number & vol_plus, ///< volume of the 'positive' part of the element
		number & vol_minus ///< volume of the 'negative' part of the element
	);
};
		
/**
 * Specialization of LSElementSize for prisms
 * \see LSElementSize
 */
template <int WDim>
class LSElementSize<ReferencePrism, WDim>
{
public:
	
	typedef ReferencePrism ref_element_type; ///< type of the reference element
	static const int dim = WDim; ///< dimension of the world
	
public:

///	computation of the volumes (for a prism)
	/**
	 * Return value: positive or negative if the entire element is in the corresponding
	 * subdomain; zero otherwise.
	 */
	static int compute
	(
		const MathVector<WDim> * corner, ///< coordinates of the corners
		const number * lsf, ///< values of the LSF at the corners
		number & vol_plus, ///< volume of the 'positive' part of the element
		number & vol_minus ///< volume of the 'negative' part of the element
	);
};

} // end namespace LevelSet
} // end namespace ug

// include implementation
#include "ls_volume_impl.h"

#endif // __H__UG__PLUGINS__LEVEL_SET__LS_VOLUME_H__

/* End of File */
