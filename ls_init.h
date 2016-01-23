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
 * Initialization of the level-set function.
 */

#ifndef __H__UG__PLUGINS__LEVEL_SET__LS_INIT_H__
#define __H__UG__PLUGINS__LEVEL_SET__LS_INIT_H__

#include <vector>

// ug4 headers
#include "common/common.h"
#include "lib_grid/algorithms/heightfield_util.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug{
namespace LevelSet{

/**
 * Initialization of a level-set function from a height raster.
 * 
 * The raster specifies the levels (heights, negated depths) of the level-set
 * function possible relatively to a specified top subset.
 *
 * \tparam	TGridFunc	type of the grid function
 */
template <typename TGridFunc>
class LSFbyRaster
{
///	grid function type
	typedef TGridFunc grid_func_type;
	
///	type of the domain
	typedef typename TGridFunc::domain_type domain_type;
	
///	type of the algebra
	typedef typename TGridFunc::algebra_type algebra_type;
	
///	grid type
	typedef typename domain_type::grid_type grid_type;
	
/// type of the position accessor
	typedef typename domain_type::position_attachment_type position_attachment_type;
	
/// type of the position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	
///	world dimension
	static const int dim = TGridFunc::dim;
	
/// vertex base iterator
	typedef typename TGridFunc::template traits<Vertex>::const_iterator VertexConstIterator;
	
	typedef lg_ntree<dim, dim, Triangle> top_tracer_tree_t;
	
	typedef RayElemIntersectionRecord<Triangle*> top_intersection_record_t;
	
public:
	
///	class constructor
	LSFbyRaster
	(
		const char * raster_file ///< ASCII raster file
	)
	:	m_bRelative (false), m_rt_tol (1e-6), m_bDefaultTop (false)
	{
		LoadHeightfieldFromASC (m_hfRaster, raster_file);
	}
	
///	computes the level-set function
	void interpolate_to
	(
		SmartPtr<TGridFunc> spLSF ///< the level-set function to compute
	);
	
///	sets the top subset
	void set_relative_to
	(
		const char * top_ss_names ///< names of the top subsets
	)
	{
		m_top_ss_names = top_ss_names;
		m_bRelative = true;
	}
	
///	sets the tolerance for getting the position of the top (only used if m_bRelative)
	void set_rel_tolerance
	(
		number tolerance ///< tolerance of the ray tracing
	)
	{
		m_rt_tol = tolerance;
	}
	
///	sets default position of the top (only used if m_bRelative)
	void set_rel_default
	(
		number default_top_z ///< default value if the ray tracing fails
	)
	{
		m_rt_default = default_top_z;
		m_bDefaultTop = true;
	}
	
private:

///	creates and fills the top tracer tree
	SmartPtr<top_tracer_tree_t> create_tt_tree
	(
		domain_type & domain ///< the domain
	);
	
/// Gets the minimum z-coordinate of the top subset over a given point
	number get_min_top_z
	(
		const MathVector<dim> & over, ///< to look over this point
		const top_tracer_tree_t & tt_tree, ///< the top tracer tree
		std::vector<top_intersection_record_t> & top_intersection_records ///< array to store all the intersections
	);
	
private:

///	the raster of the heights
	Heightfield m_hfRaster;
	
///	whether relative to the top
	bool m_bRelative;
	
///	top subsets (if specified)
	std::string m_top_ss_names;
	
///	tolerance for the ray tracing (to the top surfaces)
	number m_rt_tol;
	
///	default position of the top (if the ray tracing fails) - used only if m_bRelative
	number m_rt_default;
	
///	whether the default position of the top is specified
	bool m_bDefaultTop;
};

} // end namespace LevelSet
} // end namespace ug

// include implementation
#include "ls_init_impl.h"

#endif // __H__UG__PLUGINS__LEVEL_SET__LS_INIT_H__

/* End of File */
