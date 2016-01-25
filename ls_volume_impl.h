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
 * Measuring volumes of the parts separaterd by the level set: implementation.
 */

#include "lib_disc/reference_element/reference_element_traits.h"
#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif

namespace ug{
namespace LevelSet{

/*---- Class 'LSVolume': ----*/

/**
 * Initializes the subsets group to restrict the computation to the specified
 * subsets.
 */
template <typename TGridFunc>
void LSVolume<TGridFunc>::on_subsets
(
	const char * ss_names ///< names of the subsets
)
{
	m_ssGrp.set_subset_handler (m_spLSF->domain()->subset_handler ());
	m_ssGrp.add (TokenizeString (ss_names));
}

/**
 * Sums up the (partial) volumes of all the elements of all the types.
 */
template <typename TGridFunc>
void LSVolume<TGridFunc>::compute ()
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	sum up all the volumes
	m_volume_plus = m_volume_minus = 0;
	boost::mpl::for_each<ElemList> (AddVolumes (this));
	
#ifdef UG_PARALLEL
//	sum up the volumes from different processes
	pcl::ProcessCommunicator procComm;
	m_volume_minus = procComm.allreduce (m_volume_minus, PCL_RO_SUM);
	m_volume_plus = procComm.allreduce (m_volume_plus, PCL_RO_SUM);
#endif
}

/**
 * Sums up the (partial) volumes of all the elements of one type.
 */
template <typename TGridFunc>
template <typename TElem>
void LSVolume<TGridFunc>::add_volumes_of_all ()
{
	typedef typename grid_func_type::template traits<TElem>::const_iterator ElemIter;
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_t;

	static const size_t num_corners = ref_elem_t::numCorners;
	
	const grid_func_type & lsf = * m_spLSF;
	const position_accessor_type & aaPos = lsf.domain()->position_accessor ();
	std::vector<DoFIndex> ind (1);
	MathVector<dim> corners [num_corners];
	number lsf_values [num_corners];
	
	for (int si = 0; si < lsf.num_subsets (); si++)
	{
		if (m_ssGrp.subset_handler().valid () && ! m_ssGrp.contains (si))
			continue; // skip this subset: it is not mentioned in the specified list
		
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
					UG_THROW ("LSVolume: Not a scalar grid function for the LSF!");
				lsf_values [i] = DoFRef (lsf, ind [0]);
			}
		
		//	compute the volumes
			number vol_plus, vol_minus;
			LSElementSize<ref_elem_t, dim>::compute (corners, lsf_values, vol_plus, vol_minus);
			m_volume_plus += vol_plus; m_volume_minus += vol_minus;
			
			/*-- For debugging only: --*
			number test_vol_plus, test_vol_minus;
			number test_vol_max = (vol_plus > vol_minus)? vol_plus : vol_minus;
			for (size_t i = 0; i < num_corners; i++)
				lsf_values [i] = - lsf_values [i];
			LSElementSize<ref_elem_t, dim>::compute (corners, lsf_values, test_vol_minus, test_vol_plus);
			if (std::abs (test_vol_minus - vol_minus) / test_vol_max >= 1e-8)
			{
				UG_LOG ("---- LSVolume: Inconsistent values of the volumes:\n");
				UG_LOG ("-- V_(-) = " << vol_minus << " or " << test_vol_minus << ", diff = " << vol_minus - test_vol_minus << "\n");
				for (size_t i = 0; i < num_corners; i++)
				{
					UG_LOG ("-- co " << i << ": " << corners[i] << ", lsf = " << - lsf_values [i] << "\n");
				}
				UG_THROW ("LSVolume: Inaccurate computation of the volumes.");
			}
			 *--*/
		}
	}
}

/*---- Class 'LSElementSize': ----*/

/**
 * Specialization of 'LSElementSize<>::compute' for edges
 */
template <int WDim>
void LSElementSize<ReferenceEdge, WDim>::compute
(
	const MathVector<WDim> * corner, ///< coordinates of the corners
	const number * lsf, ///< values of the LSF at the corners
	number & vol_plus, ///< volume of the 'positive' part of the element
	number & vol_minus ///< volume of the 'negative' part of the element
)
{
	number vol_0, vol_1;
	
	if (lsf[0] * lsf[1] > 0)
	{
		vol_0 = ElementSize<ref_element_type, dim> (corner);
		vol_1 = 0;
	}
	else
	{
		number t = lsf[1] / (lsf[1] - lsf[0]);
		MathVector<dim> pnt[3];
		pnt[0] = corner[0];
		VecScaleAdd (pnt[1], t, corner[0], 1 - t, corner[1]);
		pnt[2] = corner[1];
		vol_0 = ElementSize<ref_element_type, dim> (pnt);
		vol_1 = ElementSize<ref_element_type, dim> (pnt + 1);
	}
	
	if (lsf[0] >= 0)
	{
		vol_plus = vol_0; vol_minus = vol_1;
	}
	else
	{
		vol_minus = vol_0; vol_plus = vol_1;
	}
}

/**
 * Specialization of 'LSElementSize<>::compute' for triangles
 */
template <int WDim>
void LSElementSize<ReferenceTriangle, WDim>::compute
(
	const MathVector<WDim> * corner, ///< coordinates of the corners
	const number * lsf, ///< values of the LSF at the corners
	number & vol_plus, ///< volume of the 'positive' part of the element
	number & vol_minus ///< volume of the 'negative' part of the element
)
{
	number vol = ElementSize<ref_element_type, dim> (corner);
	
//	look for 'positive' and 'negative' corners
	int i_pos = -1, i_neg = -1;
	size_t num_neg = 0;
	for (int i = 0; i < 3; i++)
		if (lsf[i] < 0)
		{
			i_neg = i;
			num_neg++;
		}
		else
			i_pos = i;
		
//	check if all are positive or all are negative
	if (num_neg == 0)
	{
		vol_plus = vol; vol_minus = 0;
		return;
	}
	if (num_neg == 3)
	{
		vol_plus = 0; vol_minus = vol;
		return;
	}
	
//	look for the cut corner
	int i_0 = (num_neg == 1)? i_neg : i_pos;
	int i_1 = (i_0 + 1) % 3, i_2 = (i_0 + 2) % 3;

//	compose a new triangle
	number t;
	MathVector<dim> cut_tri [3];
	cut_tri [0] = corner [i_0];
	t = lsf[i_1] / (lsf[i_1] - lsf[i_0]);
	VecScaleAdd (cut_tri [1], t, corner[i_0], 1 - t, corner[i_1]);
	t = lsf[i_2] / (lsf[i_2] - lsf[i_0]);
	VecScaleAdd (cut_tri [2], t, corner[i_0], 1 - t, corner[i_2]);
	
//	get the cut volume
	number cut_vol = ElementSize<ReferenceTriangle, dim> (cut_tri);
	
//	get the volumes
	if (i_0 == i_neg) // i_0 is either i_pos or i_neg
	{
		vol_minus = cut_vol; vol_plus = vol - cut_vol;
	}
	else
	{
		vol_plus = cut_vol; vol_minus = vol - cut_vol;
	}
}

/**
 * Specialization of 'LSElementSize<>::compute' for tetrahedra
 */
template <int WDim>
void LSElementSize<ReferenceTetrahedron, WDim>::compute
(
	const MathVector<WDim> * corner, ///< coordinates of the corners
	const number * lsf, ///< values of the LSF at the corners
	number & vol_plus, ///< volume of the 'positive' part of the element
	number & vol_minus ///< volume of the 'negative' part of the element
)
{
	number vol = ElementSize<ref_element_type, dim> (corner);
	
//	look for 'positive' and 'negative' corners
	int i_pos[2], i_neg[2];
	i_pos[0] = -1; i_pos[1] = -1;
	i_neg[0] = -1; i_neg[1] = -1;
	size_t num_neg = 0;
	for (int i = 0; i < 4; i++)
		if (lsf[i] < 0)
		{
			if (i_neg[0] < 0) i_neg[0] = i; else i_neg[1] = i;
			num_neg++;
		}
		else
		{
			if (i_pos[0] < 0) i_pos[0] = i; else i_pos[1] = i;
		}
	
//	check if all are positive or all are negative
	if (num_neg == 0)
	{
		vol_plus = vol; vol_minus = 0;
		return;
	}
	if (num_neg == 4)
	{
		vol_plus = 0; vol_minus = vol;
		return;
	}
	
//	if only one corner is cut out
	if (num_neg == 1 || num_neg == 3)
	{
		int i_0 = (num_neg == 1)? i_neg[0] : i_pos[0];
		MathVector<dim> cut_tet [4];
		for (int i = 0; i < 4; i++)
			if (i == i_0)
				cut_tet [i] = corner [i];
			else
			{
				number t = lsf[i] / (lsf[i] - lsf[i_0]);
				VecScaleAdd (cut_tet [i], t, corner[i_0], 1 - t, corner[i]);
			}
		// Note: The orientation of cut_tet should be the same as for the original shape, i.e. correct.
		number cut_vol = ElementSize<ReferenceTetrahedron, dim> (cut_tet);
		if (i_0 == i_neg[0]) // i_0 is either i_pos[0] or i_neg[0]
		{
			vol_minus = cut_vol; vol_plus = vol - cut_vol;
		}
		else
		{
			vol_plus = cut_vol; vol_minus = vol - cut_vol;
		}
		
		return;
	}
	
//	if two corners at every side: two prisms
	// We make the triangle containing i_neg[0] to be the botton.
	// To this end, we order i_pos[0] and i_pos[1] in such a way that
	// the new prism is properly oriented.
	int side_idx = tet_rules::FACE_FROM_VRTS [i_neg[0]] [i_pos[0]] [i_pos[1]];
	const int * side_co = tet_rules::FACE_VRT_INDS [side_idx];
	int k;
	for (k = 0; k < 3; k++)
		if (side_co [k] == i_neg[0]) break;
	UG_ASSERT (k < 3, "LSElementSize<ReferenceTetrahedron, WDim>::compute: internal error.");
	i_pos[0] = side_co [(k + 1) % 3];
	i_pos[1] = side_co [(k + 2) % 3];
	
	// Now compute the coordinates of the corners of the prism
	MathVector<dim> neg_prism [6];
	for (int k = 0; k < 2; k++)
	{
		int i_0 = i_neg [k];
		neg_prism [3 * k] = corner [i_0];
		for (int l = 0; l < 2; l++)
		{
			int i = i_pos [l];
			number t = lsf[i] / (lsf[i] - lsf[i_0]);
			VecScaleAdd (neg_prism [3 * k + l + 1], t, corner[i_0], 1 - t, corner[i]);
		}
	}
	vol_minus = ElementSize<ReferencePrism, dim> (neg_prism);
	vol_plus = vol - vol_minus;
}

/**
 * Specialization of 'LSElementSize<>::compute' for prisms
 */
template <int WDim>
void LSElementSize<ReferencePrism, WDim>::compute
(
	const MathVector<WDim> * corner, ///< coordinates of the corners
	const number * lsf, ///< values of the LSF at the corners
	number & vol_plus, ///< volume of the 'positive' part of the element
	number & vol_minus ///< volume of the 'negative' part of the element
)
{
	number vol = ElementSize<ref_element_type, dim> (corner);
	
//	look for 'positive' and 'negative' corners
	int i_pos[3], i_neg[3];
	i_pos[0] = -1; i_pos[1] = -1; i_pos[2] = -1;
	i_neg[0] = -1; i_neg[1] = -1; i_neg[2] = -1;
	size_t num_neg = 0;
	for (int i = 0; i < 6; i++)
		if (lsf[i] < 0) // note that we consider 0 as a positive value! This is used below!
		{
			if (i_neg[0] < 0) i_neg[0] = i;
			else if (i_neg[1] < 0) i_neg[1] = i;
			else i_neg[2] = i;
			num_neg++;
		}
		else
		{
			if (i_pos[0] < 0) i_pos[0] = i;
			else if (i_pos[1] < 0) i_pos[1] = i;
			else i_pos[2] = i;
		}
	/*
	 * Note that i_neg[] and i_pos[] are now ordered in ascending order.
	 */
	
//	check if all are positive or all are negative
	if (num_neg == 0)
	{
		vol_plus = vol; vol_minus = 0;
		return;
	}
	if (num_neg == 6)
	{
		vol_plus = 0; vol_minus = vol;
		return;
	}
	
//	if only one corner is cut out, i.e. a tetrahedron is cut out
	if (num_neg == 1 || num_neg == 5)
	{
		int i_0 = (num_neg == 1)? i_neg[0] : i_pos[0];
		int shift = (i_0 < 3)? 0 : 3;
		MathVector<dim> cut_tet [4];
		for (int k = 0; k < 3; k++) // loop the corners of the base where i_0 lies
		{
			int i = k + shift;
			if (i == i_0)
				cut_tet [k] = corner [i];
			else
			{
				number t = lsf[i] / (lsf[i] - lsf[i_0]);
				VecScaleAdd (cut_tet [k], t, corner[i_0], 1 - t, corner[i]);
			}
		}
		// For the base of the tetrahedron, we used the same orientation
		// as for the base of the prism. This is OK if that was the bottom
		// base of the prism. For the top base, we should invert the orientation:
		if (shift != 0)
		{
			cut_tet[3] = cut_tet[2]; // cut_tet[3] is here a temporary variable
			cut_tet[2] = cut_tet[1];
			cut_tet[1] = cut_tet[3];
		}
		int i_1 = (i_0 + 3) % 6; // connected corner of the opposite base
		number t_1 = lsf[i_1] / (lsf[i_1] - lsf[i_0]);
		VecScaleAdd (cut_tet [3], t_1, corner[i_0], 1 - t_1, corner[i_1]);
		number cut_vol = ElementSize<ReferenceTetrahedron, dim> (cut_tet);
		if (i_0 == i_neg[0]) // i_0 is either i_pos[0] or i_neg[0]
		{
			vol_minus = cut_vol; vol_plus = vol - cut_vol;
		}
		else
		{
			vol_plus = cut_vol; vol_minus = vol - cut_vol;
		}
		
		return;
	}
	
//	if two corners at every side: 2 tetrahedra or a prism cut out
	if (num_neg == 2 || num_neg == 4)
	{
		int i_0, i_1;
		if (num_neg == 2)
		{
			i_0 = i_neg[0]; i_1 = i_neg[1];
		}
		else
		{
			i_0 = i_pos[0]; i_1 = i_pos[1];
		}
		/*  Here, i_0 < i_1 (s. the remark above). */
		
	// Case 1: 2 corners of one base are cut out; this means, a prism is cut out
		if (i_1 < 3 || i_0 >= 3)
		{
			int shift = (i_0 < 3)? 0 : 3;
			int i_2 = -1; // the 3rd corner of the base
			for (int k = 0; k < 3; k++)
			{
				int i_2 = k + shift;
				if (i_2 != i_0 && i_2 != i_1) break;
			}
			int i_0a = (i_0 + 3) % 6; // corner connected to i_0 on the opposite base
			int i_1a = (i_1 + 3) % 6; // corner connected to i_1 on the opposite base
			
			MathVector<dim> cut_prism [6];
			number t;
			
			// Remark: We construct the prism in such a way that the 0th
			// side (with corners 0, 1, 4, 3) is planar (to be consistent
			// with the computation of the volume in geometry_util.h).
			// This can be achieved by taking this side from the triangular
			// base of the original prism.
			// Note however that we should provide the correct orientation
			// of the prism. Thus we should provide (i_0 + 1) % 3 == i_1
			// if the bottom base of the original prism is cut, and
			// i_0 == (i_1 + 1) % 3 for the top base. To this end, we swap
			// i_0 with i_1 (and i_0a with i_1a respectively) if necessary.
			if ((i_1 < 3 && (i_0 + 1) % 3 != i_1) || (i_0 > 3 && i_0 != (i_1 + 1) % 3))
			{
				int i;
				i = i_0; i_0 = i_1; i_1 = i;
				i = i_0a; i_0a = i_1a; i_1a = i;
			}
			
			cut_prism [0] = corner[i_0];
			t = lsf[i_2] / (lsf[i_2] - lsf[i_0]);
			VecScaleAdd (cut_prism [1], t, corner[i_0], 1 - t, corner[i_2]);
			t = lsf[i_0a] / (lsf[i_0a] - lsf[i_0]);
			VecScaleAdd (cut_prism [2], t, corner[i_0], 1 - t, corner[i_0a]);
			
			cut_prism [3] = corner[i_1];
			t = lsf[i_2] / (lsf[i_2] - lsf[i_1]);
			VecScaleAdd (cut_prism [4], t, corner[i_1], 1 - t, corner[i_2]);
			t = lsf[i_1a] / (lsf[i_1a] - lsf[i_1]);
			VecScaleAdd (cut_prism [5], t, corner[i_1], 1 - t, corner[i_1a]);
			
			number cut_vol = ElementSize<ReferencePrism, dim> (cut_prism);
			if (num_neg == 2)
			{
				vol_minus = cut_vol;
				vol_plus = vol - cut_vol;
			}
			else
			{
				vol_plus = cut_vol;
				vol_minus = vol - cut_vol;
			}
			
			return;
		}
		
	// Case 2: one corner of each base are cut out, and they are connected; a prism is cut out
		if (i_1 - i_0 == 3)
		{
			// We construct a new (cut out) prism and try to keep the 0th
			// side as planar as possible: see the previous remark. For this,
			// we take a part of a side of the original prism as a side of
			// the new one.
			MathVector<dim> cut_prism [6];
			cut_prism [0] = corner [i_0];
			cut_prism [3] = corner [i_1];
			{
				int i = (i_0 + 1) % 3; int j = i + 3;
				number t;
				t = lsf[i] / (lsf[i] - lsf[i_0]);
				VecScaleAdd (cut_prism [1], t, corner[i_0], 1 - t, corner[i]);
				t = lsf[j] / (lsf[j] - lsf[i_1]);
				VecScaleAdd (cut_prism [4], t, corner[i_1], 1 - t, corner[j]);
			}
			{
				int i = (i_0 + 2) % 3; int j = i + 3;
				number t;
				t = lsf[i] / (lsf[i] - lsf[i_0]);
				VecScaleAdd (cut_prism [2], t, corner[i_0], 1 - t, corner[i]);
				t = lsf[j] / (lsf[j] - lsf[i_1]);
				VecScaleAdd (cut_prism [5], t, corner[i_1], 1 - t, corner[j]);
			}
			number cut_vol = ElementSize<ReferencePrism, dim> (cut_prism);
			if (num_neg == 2)
			{
				vol_minus = cut_vol;
				vol_plus = vol - cut_vol;
			}
			else
			{
				vol_plus = cut_vol;
				vol_minus = vol - cut_vol;
			}
			
			return;
		}
		
	//	Case 3: The corners belong to different bases but are not connected.
	//	We consider two separated tetrahedra. They do not hinder each other.
	//	We call the function recursively.
		number tmp_lsf [6], tmp_vol;
		memcpy (tmp_lsf, lsf, 6 * sizeof (number));
		if (num_neg == 2)
		{
		//	2 negative corners
			number neg_vol;
			
			tmp_lsf [i_neg[0]] = 1; // correct the first negative corner
			LSElementSize<ReferencePrism, dim>::compute (corner, tmp_lsf, tmp_vol, neg_vol);
			vol_minus = neg_vol;
			
			tmp_lsf [i_neg[0]] = lsf [i_neg[0]];
			tmp_lsf [i_neg[1]] = 1; // correct the second negative corner
			LSElementSize<ReferencePrism, dim>::compute (corner, tmp_lsf, tmp_vol, neg_vol);
			vol_minus += neg_vol;
			
			vol_plus = vol - vol_minus;
			return;
		}
		if (num_neg == 4)
		{
		//	2 positive corners
			number pos_vol;
			
			tmp_lsf [i_pos[0]] = -1; // correct the first positive corner
			LSElementSize<ReferencePrism, dim>::compute (corner, tmp_lsf, pos_vol, tmp_vol);
			vol_plus = pos_vol;
			
			tmp_lsf [i_pos[0]] = lsf [i_pos[0]];
			tmp_lsf [i_pos[1]] = 1; // correct the second positive corner
			LSElementSize<ReferencePrism, dim>::compute (corner, tmp_lsf, pos_vol, tmp_vol);
			vol_plus += pos_vol;
			
			vol_minus = vol - vol_plus;
			return;
		}
	}
	
//	if three corners at every side: two prisms or "a very complicated case"
	if ((i_neg[0] == 0 && i_neg[2] == 2) || (i_neg[0] == 3 && i_neg[2] == 5))
	{
	//	Case 1: 2 prisms: every base has its own sign
		number cut_vol;
		MathVector<dim> cut_prism [6];
		for (int i = 0; i < 3; i++) cut_prism [i] = corner [i];
		for (int i = 0; i < 3; i++)
		{
			int j = i + 3;
			number t = lsf[j] / (lsf[j] - lsf[i]);
			VecScaleAdd (cut_prism [j], t, corner [i], 1 - t, corner [j]);
		}
		cut_vol = ElementSize<ReferencePrism, dim> (cut_prism);
		
		if (i_neg[0] == 0)
		{
			vol_minus = cut_vol;
			vol_plus = vol - cut_vol;
		}
		else
		{
			vol_plus = cut_vol;
			vol_minus = vol - cut_vol;
		}
		
		return;
	}
	
	// Case 2: The "very complicated case": Two corners of one base and one
	// of the other one are cut.
	// Only one of the three sides connecting the bases is shared by
	// a negative and a positive corner. Denote this side by [A, B], and its
	// point, where the LSF is 0, by S. (Without loss of generality, we
	// assume that LSF[A] < 0.) We divide the whole prism into two parts:
	// a smaller prism and a tetrahedron. The smaller prism is obtained from
	// the original one by shifting the corner A into S. The tetrahedron is
	// the difference of these prisms. After that, we call the function
	// recursively to compute the volumes for the smaller prism. It has
	// now one negative corner less than the original one (becase 0 is
	// considered as a positive number here), so that this case is implemented
	// above.
	MathVector<dim> part_prism [6];
	number part_lsf [6];
	int base_i = -1;
	
	for (int i = 0; i < 3; i++)
	{
		int j = i + 3;
		part_prism [i] = corner [i]; part_lsf [i] = lsf [i];
		part_prism [j] = corner [j]; part_lsf [j] = lsf [j];
		if ((lsf[i] >= 0 && lsf[j] < 0)
				|| (lsf[i] < 0 && lsf[j] >= 0)) // note: 0 is considered as a positive value
		{
			base_i = (lsf[i] < 0)? i : j; // the corner to shift
			part_lsf[base_i] = 0; // we make it to be a positive corner!
			number t = lsf[j] / (lsf[j] - lsf[i]);
			VecScaleAdd (part_prism [base_i], t, corner [i], 1 - t, corner [j]);
		}
	}
	UG_ASSERT (base_i > 0, "LSElementSize<ReferencePrism, dim>::compute: internal error");
	number part_prism_vol_plus, part_prism_vol_minus;
	LSElementSize<ReferencePrism, dim>::compute (part_prism, part_lsf,
		part_prism_vol_plus, part_prism_vol_minus);
	
	MathVector<dim> part_tet [4];
	int shift = (base_i < 3)? 0 : 3;
	for (int k = 0; k < 3; k++)
	{
		part_tet [k] = corner [k + shift];
		part_lsf [k] = lsf [k + shift];
	}
	// As before, we must provide the correct orientation of the tetrahedron.
	// To this end, we change the current orientation, if base_i is on the
	// top base of the prism. (Cf. the case of one cut corner.)
	if (shift != 0)
	{
		part_tet[3] = part_tet[2]; part_lsf[3] = part_lsf[2];
		part_tet[2] = part_tet[1]; part_lsf[2] = part_lsf[1];
		part_tet[1] = part_tet[3]; part_lsf[1] = part_lsf[3];
	}
	part_tet [3] = part_prism [base_i];
	part_lsf [3] = 0;
	number part_tet_vol_plus, part_tet_vol_minus;
	LSElementSize<ReferenceTetrahedron, dim>::compute (part_tet, part_lsf,
		part_tet_vol_plus, part_tet_vol_minus);
	
	vol_plus = part_prism_vol_plus + part_tet_vol_plus;
	vol_minus = part_prism_vol_minus + part_tet_vol_minus;
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */
