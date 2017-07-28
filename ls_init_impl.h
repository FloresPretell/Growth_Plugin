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

//	FOR DEBUGGING ONLY
// #include "lib_grid/file_io/file_io.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/gather_grid.h"
#include "lib_grid/parallelization/distributed_grid.h"
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
	position_accessor_type aaPos = spLSF->domain()->position_accessor ();
	std::vector<DoFIndex> ind (1);

//	top tracer tree
	SmartPtr<z_ray_tracer_t> sp_top_z;
	if (m_bRelative)
	{
		sp_top_z = SmartPtr<z_ray_tracer_t>
			(
				new z_ray_tracer_t (spLSF->domain (), m_top_ss_names)
			);
		sp_top_z->init (m_rt_gl);
	}
		
//	interpolate the values
	for (VertexConstIterator iter = spLSF->template begin<Vertex> ();
							iter != spLSF->template end<Vertex> (); ++iter)
	{
		Vertex * vrt = *iter;
		
	//	interpolate the raster
		MathVector<dim> coord = aaPos [vrt];
		typename raster_t::Coordinate rc;
		for (int i = 0; i < dim - 1; ++i)
			rc[i] = coord[i];
		number raster_val = m_raster.interpolate(rc, 1);
		
	//	correct the raster if relative
		if (sp_top_z.valid ())
		{
			number top_z;
			if (sp_top_z->get_min_at (coord, m_rt_tol, top_z))
				raster_val += top_z;
			else if (m_bDefaultTop)
				raster_val += m_rt_default;
			else
				UG_THROW ("LSFbyRaster: Point " << coord << " is not covered by the top subset.");
		}
	
	//	get indices of the dofs and set the value
		spLSF->inner_dof_indices (vrt, 0, ind);
		DoFRef (*spLSF, ind[0]) = coord [dim - 1] - raster_val;
	}
}

/**
 * Creates and fills the tree with the elements of the 'top' subset.
 */
template <typename TGridFunc>
void LSFbyRaster<TGridFunc>::z_ray_tracer_t::init
(
	int grid_level
)
{
	MultiGrid & mg = * m_sp_domain->grid ();
	MGSubsetHandler & sh = * m_sp_domain->subset_handler ();
	std::vector<Face*> top_faces;
	
#ifndef UG_PARALLEL
	
//	get all the elements to collect
	for (size_t i = 0; i < m_top_ss_grp.size (); i++)
	{
		int si = m_top_ss_grp [i];
		
		if (grid_level >= 0) // if the grid level for the top is specified
			for (FaceIterator it = sh.begin<Face> (si, grid_level);
									it != sh.end<Face> (si, grid_level); ++it)
				top_faces.push_back (*it);
		else
			for (int lvl = 0; lvl < (int) sh.num_levels(); lvl++)
				for (FaceIterator it = sh.begin<Face> (si, lvl);
										it != sh.end<Face> (si, lvl); ++it)
				{
					Face * t = *it;
					if (! mg.has_children (t))
						top_faces.push_back (t);
				}
	}
	
#else
	
	DistributedGridManager* dgm = mg.distributed_grid_manager();

//	select the faces on the top
	Selector sel (mg);
	for (size_t i = 0; i < m_top_ss_grp.size (); i++)
	{
		int si = m_top_ss_grp [i];
		
		if (grid_level >= 0) // if the grid level for the top is specified
			for (FaceIterator it = sh.begin<Face> (si, grid_level);
									it != sh.end<Face> (si, grid_level); ++it)
				sel.select (*it);
		else
			for (int lvl = 0; lvl < (int) sh.num_levels(); lvl++)
				for (FaceIterator it = sh.begin<Face> (si, lvl);
										it != sh.end<Face> (si, lvl); ++it)
				{
					Face * t = *it;
					if (! (mg.has_children (t) || (dgm && dgm->is_ghost(t))))
						sel.select (t);
				}
	}
	
//	copy the top faces into a new grid
	GridDataSerializationHandler serializer;
	serializer.add
		(GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (mg, m_sp_domain->position_attachment ()));
		
	GridDataSerializationHandler deserializer;
	deserializer.add
		(GeomObjAttachmentSerializer<Vertex, position_attachment_type>::create (m_top_grid, m_sp_domain->position_attachment ()));

	AllGatherGrid (m_top_grid, sel, serializer, deserializer);

	for (FaceIterator it = m_top_grid.begin<Face> (); it != m_top_grid.end<Face> (); ++it)
		top_faces.push_back (*it);
	
	// UG_LOG("DEBUG: SAVING allgathered m_top_grid to file in ls_init...\n");
	// SaveGridToFile(m_top_grid, mkstr("top_grid_p" << pcl::ProcRank() << ".ugx").c_str(),
	//                m_sp_domain->position_attachment());
#endif // UG_PARALLEL
	
//	compose the tree
	m_top_tracer_tree.create_tree (top_faces.begin (), top_faces.end ());
}

/**
 * Gets the minimum z-coordinate of the top subset over a given point
 */
template <typename TGridFunc>
bool LSFbyRaster<TGridFunc>::z_ray_tracer_t::get_min_at
(
	const MathVector<dim> & over, ///< to look over this point
	number tolerance, ///< the tolerance of the ray tracer
	number & z ///< the result
)
{
//	find all the intersections
	MathVector<dim> up_dir;
	up_dir = 0;
	up_dir [dim - 1] = 1;
	m_top_intersection_records.clear ();
	RayElementIntersections (m_top_intersection_records, m_top_tracer_tree, over, up_dir, tolerance);
	
//	check if there are intersections at all
	if (m_top_intersection_records.size () == 0)
		return false;
	
//	find the lowest point
	MathVector<dim> x = PointOnRay (over, up_dir, m_top_intersection_records[0].smin);
	number z_min = x [dim - 1];
	for (size_t i = 1; i < m_top_intersection_records.size (); i++)
	{
		top_intersection_record_t & r = m_top_intersection_records [i];
		x = PointOnRay (over, up_dir, r.smin);
		if (x [dim - 1] < z_min)
			z_min = x [dim - 1];
	}
	
	z = z_min;
	return true;
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */