/*
 * Tools for determining the position of the level set
 *
 * Created on Jul. 10, 2015 by D. Logashenko
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET_POSITION_H__
#define __H__UG__PLUGINS__LEVEL_SET_POSITION_H__

#include <vector>

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "lib_disc/reference_element/element_list_traits.h"

namespace ug {
namespace LevelSet {

/**
 * Class for measuring the z-coordinates of the level set
 */
template <typename TGridFunction>
class LSPositionZ
{
public:

/// grid function type
	typedef TGridFunction grid_func_type;
	
///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

///	generic element type
	typedef typename TGridFunction::element_type elem_type;
	
///	world dimension
	static const int dim = TGridFunction::dim;
	
/// max. number of the corners of the elements in the domain
	static const size_t maxNumCorners = (size_t)
		element_list_traits<typename domain_traits<dim>::DimElemList>::maxCorners;
	
private:
	
///	threshold for the level-set-function
	static number lsf_threshold () {return 1e-8;}
	
///	quad- or octtree type to search the elements
	typedef lg_ntree<dim-1, dim, elem_type> tree_type;

///	intersection record type
    typedef RayElemIntersectionRecord<elem_type *> intersect_rec_type;
    
///	a structure to get the intersection points of an element
	struct elem_intersect_data
	{
		elem_type * elem; ///< the intersected element
		MathVector<dim> pnt[2]; ///< the intersection points
		
	//	constructor
		elem_intersect_data (elem_type * e, MathVector<dim> p_min, MathVector<dim> p_max)
		:	elem (e)
		{
			pnt[0] = p_min; pnt[1] = p_max;
		}
	};
    
public:

///	class constructor
	LSPositionZ
	(
		SmartPtr<grid_func_type> spLSF ///< the level-set function
	)
	:	m_spLSF (spLSF),
		m_tree (* (spLSF->domain()->grid ()), spLSF->domain()->position_attachment ())
	{}
	
///	initialize the object: find the intersected elements and prepare the tree
	void reinit ()
	{
	//	find the intersected elements
		fill_cut ();
	
	//	create the tree
		m_tree.create_tree (m_vpCut.begin (), m_vpCut.end ());
	};
	
///	get the 'z'-coordinates
	void get_height
	(
		MathVector<dim-1> xy, ///< the first coordinates (x, y)
		std::vector<number> & z, ///< to save the last coordinates (of all the intersections)
		number small_z = 1e-12 ///< a threshold for the coordinates
	)
	{
		z.clear ();
		if (! get_intersected (xy)) return;
		get_z (z, small_z);
	}
	
///	gets the largest 'z'-coordinate
	number get_height_at
	(
		const MathVector<dim> & pnt, ///< the first coordinates (x, y), die z-coordinate is neglected
		number default_value ///< returned if no intersections found
	)
	{
		std::vector<number> z;
		MathVector<dim-1> xy;
		for (size_t i = 0; i < dim-1; i++)
			xy[i] = pnt[i];
		get_height (xy, z);
		if (z.size () == 0) return default_value;
		number max_z = z[0];
		for (size_t i = 1; i < z.size (); i++)
			if (z[i] > max_z)
				max_z = z[i];
		return max_z;
	}
	
private:

///	sign of the LSF in the element
	inline int lsf_sign (size_t noc, number lsf []);
	
///	get the elements cut by the level set
	void fill_cut ();
	
///	get the elements intersected by the ray
	bool get_intersected
	(
		MathVector<dim-1> xy ///< the first coordinates of the points of the line
	);
	
///	get the point at the level set
	void get_z
	(
		std::vector<number> & z, ///< the z-coordinates of the intersections
		number small_z = 1e-12 ///< a threshold for the coordinates
	);
	
private:

///	level-set function
	SmartPtr<grid_func_type> m_spLSF;

///	the quad- or octtree:
	tree_type m_tree;
	
/// array of the intersected elements
	std::vector<elem_type *> m_vpCut;
	
/// intersection recorder
	std::vector<intersect_rec_type> m_intersectionRecords;
	
///	set of the intersections
	std::vector<elem_intersect_data> m_intersections;
};

} // namespace LevelSet
} // end namespace ug

// include implementation
#include "level_set_pos_impl.h"

#endif // __H__UG__PLUGINS__LEVEL_SET_POSITION_H__

/* End of File */
