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
 * Initialization of the level-set function: implementation
 */

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/broadcast.h" 
#endif

namespace ug{
namespace LevelSet{

/**
 * Interpolates the raster at the grid points and corrects the values by
 * the level of the top subset (if specified).
 */
template <typename TGridFunc>
void LSFbyRaster<TGridFunc>::interpolate_to
(
	SmartPtr<TGridFunc> spLSF ///< the level-set function to compute
)
{
	domain_type & domain = * spLSF->domain().get ();
	position_accessor_type aaPos = domain.position_accessor ();
	std::vector<DoFIndex> ind (1);

//	top tracer tree
	SmartPtr<top_tracer_tree_t> sp_top_tracer_tree;
	std::vector<top_intersection_record_t> top_intersection_records;
	if (m_bRelative)
		sp_top_tracer_tree = create_tt_tree (domain);
		
//	interpolate the values
	for (VertexConstIterator iter = spLSF->template begin<Vertex> ();
							iter != spLSF->template end<Vertex> (); ++iter)
	{
		Vertex * vrt = *iter;
		
	//	interpolate the raster
		MathVector<dim> coord = aaPos [vrt];
		number raster_val = m_hfRaster.interpolate (coord[0], coord[1]); //TODO: this excludes 1d!
		
	//	correct the raster if relative
		if (sp_top_tracer_tree.valid ())
			raster_val += get_min_top_z (coord, *sp_top_tracer_tree, top_intersection_records);
	
	//	get indices of the dofs and set the value
		spLSF->inner_dof_indices (vrt, 0, ind);
		DoFRef (*spLSF, ind[0]) = coord [dim - 1] - raster_val;
	}
}

/**
 * Creates and fills the tree with the elements of the 'top' subset.
 */
template <typename TGridFunc>
SmartPtr<typename LSFbyRaster<TGridFunc>::top_tracer_tree_t> LSFbyRaster<TGridFunc>::create_tt_tree
(
	domain_type & domain ///< the domain
)
{
	MultiGrid & mg = * domain.grid ();
	MGSubsetHandler & sh = * domain.subset_handler ();
	std::vector<Triangle*> tris;
	
//	parse the subset names
	SubsetGroup ssGrp (domain.subset_handler ());
	ssGrp.add (TokenizeString (m_top_ss_names));
	
#ifndef UG_PARALLEL
//	create the tree
	SmartPtr<top_tracer_tree_t> sp_top_tracer_tree (new top_tracer_tree_t (mg, domain.position_attachment ()));
	
//	get all the elements to collect
	for (size_t i = 0; i < ssGrp.size (); i++)
	{
		int si = ssGrp [i];
		for (int lvl = 0; lvl < (int) sh.num_levels(); lvl++)
		{
			for (TriangleIterator it = sh.begin<Triangle> (si, lvl);
									it != sh.end<Triangle> (si, lvl); ++it)
			{
				Triangle * t = *it;
				if (! mg.has_children (t))
					tris.push_back (t);
			}
		}
	}
#else
	const int rootProc = 0;
	
//	select the faces on the top
	Selector sel (mg);
	if (pcl::ProcRank () == rootProc)
	{
		for (size_t i = 0; i < ssGrp.size (); i++)
		{
			int si = ssGrp [i];
			for (int lvl = 0; lvl < (int) sh.num_levels(); lvl++)
			{
				for (TriangleIterator it = sh.begin<Triangle> (si, lvl);
										it != sh.end<Triangle> (si, lvl); ++it)
				{
					Triangle * t = *it;
					if (! mg.has_children (t))
						sel.select (t);
				}
			}
		}
	}
	
//	copy the top faces into a new grid
	GridDataSerializationHandler serializer;
	serializer.add
		(GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (mg, domain.position_attachment ()));
	Grid top_grid;
	BroadcastGrid (top_grid, sel, serializer, rootProc);
	for (TriangleIterator it = top_grid.begin<Triangle> (); it != top_grid.end<Triangle> (); ++it)
		tris.push_back (*it);
	
//	create the tree
	SmartPtr<top_tracer_tree_t> sp_top_tracer_tree (new top_tracer_tree_t (top_grid, domain.position_attachment ()));
#endif // UG_PARALLEL
	
//	compose the tree
	sp_top_tracer_tree->create_tree (tris.begin (), tris.end ());
	
	return sp_top_tracer_tree;
}

/**
 * Gets the minimum z-coordinate of the top subset over a given point
 */
template <typename TGridFunc>
number LSFbyRaster<TGridFunc>::get_min_top_z
(
	const MathVector<dim> & over, ///< to look over this point
	const top_tracer_tree_t & tt_tree, ///< the top tracer tree
	std::vector<top_intersection_record_t> & top_intersection_records ///< array to store all the intersections
)
{
//	find all the intersections
	MathVector<dim> up_dir;
	up_dir = 0;
	up_dir [dim - 1] = 1;
	top_intersection_records.clear ();
	RayElementIntersections (top_intersection_records, tt_tree, over, up_dir, m_rt_tol);
	
//	check if there are intersections at all
	if (top_intersection_records.size () == 0)
	{
		if (m_bDefaultTop)
			return m_rt_default;
		UG_THROW ("LSFbyRaster: Point " << over << " is not covered by the top subset.");
	}
	
//	find the lowest point
	MathVector<dim> x = PointOnRay (over, up_dir, top_intersection_records[0].smin);
	number z_min = x [dim - 1];
	for (size_t i = 1; i < top_intersection_records.size (); i++)
	{
		top_intersection_record_t & r = top_intersection_records [i];
		x = PointOnRay (over, up_dir, r.smin);
		if (x [dim - 1] < z_min)
			z_min = x [dim - 1];
	}
	
	return z_min;
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */