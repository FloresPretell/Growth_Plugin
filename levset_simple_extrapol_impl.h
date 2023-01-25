/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
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
 * Implementation of functions from levset_simple_extrapol.h
 */

// ug4 headers
#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#endif

namespace ug {
namespace LevelSet {

/*-------- class LevSetGFsimpleExtrapolation --------*/

/**
 * Gets the level-set function at corners of an element and corrects it to avoid
 * small "inside" values. The function returns a positive value if the level-set
 * function is not specified or the element is completely 'inside'. If the element
 * is intersected by the interface, the function returns 0. If the element is
 * completely outside, the function returns a negative value.
 */
template <typename TDomain, typename TAlgebra>
int LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::check_elem_lsf
(
	size_t n_co, ///< number of the corners of the element
	GridObject * pElem, ///< the element to process
	int si, ///< subset of the element
	int g_level, ///< grid level of the element
	bool use_hanging, ///< if there can be hanging nodes
	const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
	number time ///< the phisical time
)
{
	if (m_spLSF.invalid ()) return 1; // "no LSF" == "completely inside"
	
//	Hanging nodes are not supported
	if (use_hanging)
		UG_THROW ("LevSetGFsimpleExtrapolation: Hanging nodes are not supported.");
	
//	the level-set function on the level
	UG_ASSERT (g_level >= 0 && (size_t) g_level < m_vGLData.size (), "Grid level of the level-set function mismatch!");
	SmartPtr<ls_grid_func_type> spLSF = m_vGLData[g_level].lsf_on_gl;
	
//	keep the current element
	m_pElem = pElem; m_si = si; m_vCornerCoords = vCornerCoords;
	m_time = time;
	
//	get the corner values of the LSF:
	spLSF->indices (pElem, m_indLSF, false);
	m_locLSF.resize (m_indLSF);
	GetLocalVector (m_locLSF, *spLSF);

//	check if the element is inside, or outside, or intersected
	int inside, outside;
	
    inside = outside = 0;
	for (size_t co = 0; co < n_co; co++)
	{
		if (corner_inside (co))
			inside = 1;
		else
			outside = 1;
	}
	
	inside -= outside;
	
	if (inside != 0)
		return inside; // not intersected
	
//	for an intersected element, compute the center
	m_elemCenter = vCornerCoords[0];
	for (size_t co = 1; co < n_co; co++)
		m_elemCenter += vCornerCoords[co];
	m_elemCenter /= (number) n_co;
	
//	get the corners
	if (m_excl_ssg.size () != 0)
	{
		Grid::vertex_traits::secure_container corner_list;
		m_excl_ssg.subset_handler()->grid()->associated_elements (corner_list, pElem);
		if (corner_list.size () != n_co) // hanging nodes?
			UG_THROW ("LevSetGFsimpleExtrapolation: Hanging nodes are not supported. - Illegal number of vertices in an element.");
		for (size_t co = 0; co < n_co; co++)
			m_co_excluded [co] = m_excl_ssg.contains (m_excl_ssg.subset_handler()->get_subset_index (corner_list [co]));
	}
	else
		for (size_t co = 0; co < max_num_corners; co++)
			m_co_excluded [co] = false;
	
	return 0;
}

/**
 * Gets the boundary condition value at the embedded interface.
 */
template <typename TDomain, typename TAlgebra>
number LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::boundary_value
(
	size_t fct ///< the function
) const
{
	const ICData & ic_data = m_vICData[fct];
	
	if (ic_data.func.invalid ())
		return ic_data.value;
	
	number val;
	(* ic_data.func) (val, m_elemCenter, m_time, m_si);
	return val;
}

/**
 * Extrapolates all the components of the solution to the vertices behind the interface.
 *
 * This function extrapolates all the 'outer' nodal values using the value
 * at the given base corner.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::extrapolate_sol_by_lsf
(
	LocalVector& locU, ///< the solution to extrapolate
	size_t base_co ///< the base corner (to extrapolate from); note that it should be 'inside'
) const
{
	for (size_t fct = 0; fct < locU.num_all_fct (); fct++)
		if (locU.num_all_dof (fct) != TElem::NUM_VERTICES)
			UG_THROW ("LevSetGFsimpleExtrapolation:"
				" Hanging nodes are not currently supported for the ghost-fluid method");
	
	for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
	if ((! corner_inside (co)) && (! corner_excluded (co))) /* extrapolate */
	{
		number t = m_locLSF.value (0, co) / m_locLSF.value (0, base_co);
		for (size_t fct = 0; fct < locU.num_all_fct (); fct++)
			switch (m_vICData[fct].bc_kind)
			{
				case DIRICHLET_BC:
					locU.value (fct, co) = locU.value (fct, base_co) * t + boundary_value (fct) * (1 - t);
					break;
					
				case DIRICHLET_PLAIN_BC:
					locU.value (fct, co) = boundary_value (fct);
					break;
				
				case NEUMANN_0_BC:
					locU.value (fct, co) = locU.value (fct, base_co);
					break;
				
				default:
					UG_THROW ("LevSetGFsimpleExtrapolation: Unknown interface condition type.");
			}
	}
	//	else use the original values (at the corners 'inside')
}

/**
 * Extrapolates a component of the solution to the vertices behind the interface w.r.t. a base corner.
 *
 * This function extrapolates all the 'outer' nodal values of one component using the value
 * at the given base corner.
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::extrapolate_by_lsf
(
	size_t num_co, ///< number of the corners
	size_t base_co, ///< the base corner (to extrapolate from); note that it should be 'inside'
	number * u, ///< nodal values to extrapolate
	size_t fct ///< index of the function (to identify to type of the extrapolation)
) const
{
	if (fct < 0 || fct >= m_vICData.size ())
		UG_THROW ("LevSetGFsimpleExtrapolation: Wrong function index.");
	
	t_fs_bc bc_kind = m_vICData[fct].bc_kind;
	number interface_val = boundary_value (fct);
	
	for (size_t co = 0; co < num_co; co++)
	if ((! corner_inside (co)) && (! corner_excluded (co))) /* extrapolate */
		switch (bc_kind)
		{
			case DIRICHLET_BC:
				{
					number t = m_locLSF.value (0, co) / m_locLSF.value (0, base_co);
					u [co] = u [base_co] * t + interface_val * (1 - t);
				}
				break;
				
			case DIRICHLET_PLAIN_BC:
				u [co] = interface_val;
				break;
			
			case NEUMANN_0_BC:
				u [co] = u [base_co];
				break;
			
			default:
				UG_THROW ("LevSetGFsimpleExtrapolation: Unknown interface condition type.");
		}
	//	else use the original values (at the corners 'outside')
}

/**
 * Extrapolates a component of the solution to the vertices behind the interface by averaging.
 *
 * This function extrapolates every 'outer' nodal values by the average of the
 * extrapolations from the 'inner' nodes
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::extrapolate_by_lsf
(
	size_t num_co, ///< number of the corners
	number * u, ///< nodal values to extrapolate
	size_t fct ///< index of the function (to identify to type of the extrapolation)
) const
{
	if (fct < 0 || fct >= m_vICData.size ())
		UG_THROW ("LevSetGFsimpleExtrapolation: Wrong function index.");
	
	t_fs_bc bc_kind = m_vICData[fct].bc_kind;
	number interface_val = boundary_value (fct);
	
	for (size_t co = 0; co < num_co; co++)
	if ((! corner_inside (co)) && (! corner_excluded (co))) /* extrapolate */
	{
		size_t n_base_co = 0;
		u [co] = 0;
		for (size_t base_co = 0; base_co < num_co; base_co++)
		if (corner_inside (base_co)) /* use this corner to extrapolate */
		{
			switch (bc_kind)
			{
				case DIRICHLET_BC:
					{
						number t = m_locLSF.value (0, co) / m_locLSF.value (0, base_co);
						u [co] += u [base_co] * t + interface_val * (1 - t);
					}
					break;
				
				case DIRICHLET_PLAIN_BC:
					u [co] += interface_val;
					break;
			
				case NEUMANN_0_BC:
					u [co] += u [base_co];
					break;
			
				default:
					UG_THROW ("LevSetGFsimpleExtrapolation: Unknown interface condition type.");
			}
			n_base_co++;
		}
		UG_ASSERT (n_base_co != 0, "LevSetGFsimpleExtrapolation:"
			"Attempt to interpolate in an element that is not cut.");
		u [co] /= n_base_co;
	}
	//	else use the original values (at the corners 'outside')
}

/**
 * Sets the values of the vector at the non-base vertices to 0.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::clear_outer_vectors
(
	LocalVector& locD, ///< the local defect/rhs to postprocess
	size_t base_co ///< the base corner (to extrapolate from)
)
{
	for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
	if (co != base_co)
		for (size_t fct = 0; fct < locD.num_all_fct (); fct++)
			locD.value (fct, co) = 0;
}

/**
 * Eliminates the matrix connections to the vertices behind the interface.
 * After that, this function sets to 0 all the rows of the local matrix
 * except that for the base corner.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::eliminate_extrapolated
(
	LocalMatrix& locM, ///< the local matrix to postprocess
	size_t base_co ///< the base corner (to extrapolate from)
)
{
//	eliminate
	for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
	if ((! corner_inside (co)) && (! corner_excluded (co)))
	{
		number t = m_locLSF.value (0, co) / m_locLSF.value (0, base_co);
		for (size_t row_fct = 0; row_fct < locM.num_all_row_fct (); row_fct++)
			for (size_t col_fct = 0; col_fct < locM.num_all_col_fct (); col_fct++)
			{
				number & a_ij = locM.value (row_fct, base_co, col_fct, co);
				switch (m_vICData[col_fct].bc_kind)
				{
					case DIRICHLET_BC:
						locM.value (row_fct, base_co, col_fct, base_co) += a_ij * t;
						break;
					
					case DIRICHLET_PLAIN_BC:
						// noting to do: the contribution does not depend on the unknowns
						break;
				
					case NEUMANN_0_BC:
						locM.value (row_fct, base_co, col_fct, base_co) += a_ij;
						break;
				
					default:
						UG_THROW ("LevSetGFsimpleExtrapolation: Unknown interface condition type.");
				}
				a_ij = 0;
			}
	}
	
//	clear
	for (size_t row_co = 0; row_co < TElem::NUM_VERTICES; row_co++)
	if (row_co != base_co)
		for (size_t col_co = 0; col_co < TElem::NUM_VERTICES; col_co++)
			for (size_t row_fct = 0; row_fct < locM.num_all_row_fct (); row_fct++)
				for (size_t col_fct = 0; col_fct < locM.num_all_col_fct (); col_fct++)
					locM.value (row_fct, row_co, col_fct, col_co) = 0;
}

/**
 * Eliminates the matrix connections to the vertices behind the interface and
 * corrects the right-hand side. After that, this function sets to 0 all the
 * rows of the local matrix and all the entries of the local defect at the base
 * corner.
 *
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::eliminate_extrapolated
(
	LocalMatrix& locM, ///< the local matrix to postprocess
	LocalVector& locB, ///< the local right-hand side to postprocess
	size_t base_co ///< the base corner (to extrapolate from)
)
{
//	eliminate
	for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
	if ((! corner_inside (co)) && (! corner_excluded (co)))
	{
		number t = m_locLSF.value (0, co) / m_locLSF.value (0, base_co);
		for (size_t row_fct = 0; row_fct < locM.num_all_row_fct (); row_fct++)
			for (size_t col_fct = 0; col_fct < locM.num_all_col_fct (); col_fct++)
			{
				number & a_ij = locM.value (row_fct, base_co, col_fct, co);
				switch (m_vICData[col_fct].bc_kind)
				{
					case DIRICHLET_BC:
						locM.value (row_fct, base_co, col_fct, base_co) += a_ij * t;
						locB.value (row_fct, base_co) -= a_ij * boundary_value (col_fct) * (1 - t);
						break;
					
					case DIRICHLET_PLAIN_BC:
						locB.value (row_fct, base_co) -= a_ij * boundary_value (col_fct);
						break;
				
					case NEUMANN_0_BC:
						locM.value (row_fct, base_co, col_fct, base_co) += a_ij;
						break;
				
					default:
						UG_THROW ("LevSetGFsimpleExtrapolation: Unknown interface condition type.");
				}
				a_ij = 0;
			}
	}
	
//	clear
	for (size_t row_co = 0; row_co < TElem::NUM_VERTICES; row_co++)
	if (row_co != base_co)
		for (size_t row_fct = 0; row_fct < locM.num_all_row_fct (); row_fct++)
		{
			locB.value (row_fct, row_co) = 0;
			for (size_t col_co = 0; col_co < TElem::NUM_VERTICES; col_co++)
				for (size_t col_fct = 0; col_fct < locM.num_all_col_fct (); col_fct++)
					locM.value (row_fct, row_co, col_fct, col_co) = 0;
		}
}

/**
 * Eliminates the matrix connections to the vertices behind the interface and
 * corrects the right-hand side. After that, this function sets to 0 all the
 * rows of the local matrix and all the entries of the local defect at the outer
 * corners.
 *
 * REMARK: This function should be used for linear problems only!
 *
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::eliminate_extrapolated
(
	LocalMatrix& locM, ///< the local matrix to postprocess
	LocalVector& locB ///< the local right-hand side to postprocess
)
{
//	eliminate
	for (size_t base_co = 0; base_co < TElem::NUM_VERTICES; base_co++)
	if (corner_inside (base_co))
		for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
		if ((! corner_inside (co)) && (! corner_excluded (co)))
		{
			number t = m_locLSF.value (0, co) / m_locLSF.value (0, base_co);
			for (size_t row_fct = 0; row_fct < locM.num_all_row_fct (); row_fct++)
				for (size_t col_fct = 0; col_fct < locM.num_all_col_fct (); col_fct++)
				{
					number & a_ij = locM.value (row_fct, base_co, col_fct, co);
					switch (m_vICData[col_fct].bc_kind)
					{
						case DIRICHLET_BC:
							locM.value (row_fct, base_co, col_fct, base_co) += a_ij * t;
							locB.value (row_fct, base_co) -= a_ij * boundary_value (col_fct) * (1 - t);
							break;
					
						case DIRICHLET_PLAIN_BC:
							locB.value (row_fct, base_co) -= a_ij * boundary_value (col_fct);
							break;
				
						case NEUMANN_0_BC:
							locM.value (row_fct, base_co, col_fct, base_co) += a_ij;
							break;
				
						default:
							UG_THROW ("LevSetGFsimpleExtrapolation: Unknown interface condition type.");
					}
					a_ij = 0;
				}
		}
	
//	clear
	for (size_t row_co = 0; row_co < TElem::NUM_VERTICES; row_co++)
	if (! corner_inside (row_co))
		for (size_t row_fct = 0; row_fct < locM.num_all_row_fct (); row_fct++)
		{
			locB.value (row_fct, row_co) = 0;
			for (size_t col_co = 0; col_co < TElem::NUM_VERTICES; col_co++)
				for (size_t col_fct = 0; col_fct < locM.num_all_col_fct (); col_fct++)
					locM.value (row_fct, row_co, col_fct, col_co) = 0;
		}
}

/**
 * Prepares the grid functions (as data structures - not values) and the
 * restriction operators (canonical injections) on grid levels.
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::prepare_grid_levels ()
{
	m_vGLData.resize (0); // to deallocate the old data
	if (! m_spLSF.valid ()) return; // no level-set function
	
	int finest_lev = m_spLSF->grid_level().level ();
	GridLevel::ViewType view_type = m_spLSF->grid_level().type (); //TODO: Is it correct, to preserve the view type?
	SmartPtr<ApproximationSpace<domain_type> > approx_space = m_spLSF->approx_space ();
	
	if (finest_lev == GridLevel::TOP)
		finest_lev = approx_space->num_levels () - 1;
	UG_ASSERT (finest_lev >= 0, "Wrong finest grid level!");
	
//	the finest grid level (as specified by the original LSF)
	m_vGLData.resize (finest_lev + 1);
	m_vGLData[finest_lev].lsf_on_gl = m_spLSF;
	m_vGLData[finest_lev].inject = SPNULL;
	
//	coarser grid levels
	for (int fine_lev = finest_lev; fine_lev > 0; fine_lev--)
	{
		int coarse_lev = fine_lev - 1;
		GridLevel fine_gl (fine_lev, view_type);
		GridLevel coarse_gl (coarse_lev, view_type);
		
		m_vGLData[coarse_lev].lsf_on_gl = SmartPtr<ls_grid_func_type> (new ls_grid_func_type (approx_space, coarse_gl, false));
		
		m_vGLData[coarse_lev].inject = SmartPtr<projection_type> (new projection_type (approx_space));
		m_vGLData[coarse_lev].inject->set_levels (coarse_gl, fine_gl);
		m_vGLData[coarse_lev].inject->init ();
	}
}

/**
 * Projects the level-set function to the coarser grid levels
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::project_LSF ()
{
	if (m_vGLData.size () < 2) return; // nothing to project
	for (int coarse_lev = m_vGLData.size () - 2; coarse_lev >= 0; coarse_lev--)
		m_vGLData[coarse_lev].inject->do_restrict
			(* m_vGLData[coarse_lev].lsf_on_gl, * m_vGLData[coarse_lev + 1].lsf_on_gl);
}

/**
 * Sets all vertex dofs at outer vertices to 0.
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::clear_outer_values
(
	vector_type & d, ///< the vector where to set
	const DoFDistribution * dd ///< dof distribution of the grid function to reset
) const
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator t_vert_iterator;
	
//	If no LSF given, do nothing
	if (m_spLSF.invalid ()) return;
	
//	Grid level of the dof distribution, and the correct level-set function
	int level = dd->grid_level().level ();
	if (level == GridLevel::TOP)
		level = dd->multi_grid()->top_level ();
	if (level < 0 || (size_t) level >= m_vGLData.size ())
		UG_THROW ("Attempt to assemble on a grid level where the LSF is undefined.");
	SmartPtr<ls_grid_func_type> spLSF = m_vGLData[level].lsf_on_gl;

//	Arrays for the indices in the grid functions:
	std::vector<size_t> vLSFVertInd (1);
	std::vector<size_t> vDefVertInd;
	
//	Loop the vertices
	t_vert_iterator iter = dd->template begin<Vertex> ();
	t_vert_iterator iterEnd = dd->template end<Vertex> ();
	for (; iter != iterEnd; iter++)
	{
		Vertex * pVertex = *iter;
		
	//	Check if the vertex is excluded:
		if (m_excl_ssg.size () != 0 && m_excl_ssg.contains (m_excl_ssg.subset_handler()->get_subset_index (pVertex)))
			continue;
		
	//	Get the multiindex of the LSF and check the value of the LSF
		if (spLSF->inner_algebra_indices (pVertex, vLSFVertInd) != 1)
			UG_THROW ("LevSetGFsimpleExtrapolation: Non-scalar Level-Set Function.");
		if (lsf_inside (BlockRef ((* spLSF) [vLSFVertInd[0]], 0)))
			continue;
		
	//	Get the multiindices of the grid function and set the values:
		size_t n_dofs = dd->inner_algebra_indices (pVertex, vDefVertInd);
		for (size_t dof = 0; dof < n_dofs; dof++)
			d [vDefVertInd[dof]] = 0;
	}
}

/**
 * Sets all vertex dofs at outer vertices to given values.
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::set_outer_values
(
	vector_type & u, ///< the vector where to set
	const DoFDistribution * dd, ///< dof distribution of the grid function to reset
	number time ///< the physical time
)
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator t_vert_iterator;
	
//	If no LSF given, do nothing
	if (m_spLSF.invalid ()) return;

//	Grid level of the dof distribution, and the correct level-set function
	int level = dd->grid_level().level ();
	if (level == GridLevel::TOP)
		level = dd->multi_grid()->top_level ();
	if (level < 0 || (size_t) level >= m_vGLData.size ())
		UG_THROW ("LevSetGFsimpleExtrapolation: Attempt to assemble on a grid level where the LSF is undefined.");
	SmartPtr<ls_grid_func_type> spLSF = m_vGLData[level].lsf_on_gl;

	std::vector<size_t> vLSFVertInd (1);
	std::vector<DoFIndex> multInd (1);
	
	ANumber aBC;
	AUInt aNumElem;
	grid_type & grid = * (grid_type *) (dd->multi_grid().get ()); // we cancel the 'const' specifier here!
	grid.attach_to_vertices (aBC);
	grid.attach_to_vertices (aNumElem);
	Grid::VertexAttachmentAccessor<ANumber> aaBC (grid, aBC);
	Grid::VertexAttachmentAccessor<AUInt> aaNumElem (grid, aNumElem);
	typedef typename domain_traits<dim>::DimElemList AssembleElemList;
	
//	Initialize the outer vertices near the interface: Sum up the extrapolated values
	for (size_t fct = 0; fct < dd->num_fct (); fct++)
	{
	//	Prepare the attachments
		t_vert_iterator iter = dd->template begin<Vertex> ();
		t_vert_iterator iterEnd = dd->template end<Vertex> ();
		for (; iter != iterEnd; iter++)
		{
			Vertex * pVertex = *iter;
			aaBC [pVertex] = 0;
			aaNumElem [pVertex] = 0;
		}
		
	//	Sup up and count the values
		boost::mpl::for_each<AssembleElemList> (SumUpNearIfOuterValues (this, u, fct, dd, time, aaBC, aaNumElem));
#		ifdef UG_PARALLEL
		AttachmentAllReduce<Vertex> (grid, aBC, PCL_RO_SUM);
		AttachmentAllReduce<Vertex> (grid, aNumElem, PCL_RO_SUM);
#		endif
	
	//	Average the extrapolated values
		for (int si = 0; si < dd->num_subsets (); si++)
		{
		//	Check if the function index is defined in this subset
			if (! dd->is_def_in_subset (fct, si))
				continue;
			
		//	Check if the subset is excluded:
			if (m_excl_ssg.size () != 0 && m_excl_ssg.contains (si))
				continue;
	
		//	Loop the vertices
			t_vert_iterator iter = dd->template begin<Vertex> (si);
			t_vert_iterator iterEnd = dd->template end<Vertex> (si);
			for (; iter != iterEnd; iter++)
			{
				Vertex * pVertex = *iter;
			
			//	Get the multiindex of the LSF and check the value of the LSF
				if (spLSF->inner_algebra_indices (pVertex, vLSFVertInd) != 1)
					UG_THROW ("LevSetGFsimpleExtrapolation: Non-scalar Level-Set Function.");
				if (lsf_inside (BlockRef ((* spLSF) [vLSFVertInd[0]], 0)))
					continue; // we do not reset the values that are inside
					
				if (dd->inner_dof_indices (pVertex, fct, multInd) != 1)
					UG_THROW ("LevSetGFsimpleExtrapolation: More than one DoF per vertex for a component. Not the Lagrange element?");
				
				const uint nElem = aaNumElem [pVertex];
				if (nElem != 0) // if 0, the values have not been extrapolated at this vertex
				{
					const number bcVal = aaBC [pVertex];
					DoFRef (u, multInd[0]) = bcVal / nElem;
				}
				else
					DoFRef (u, multInd[0]) = 0;
			}
		}
	}
	
	grid.detach_from_vertices (aNumElem);
	grid.detach_from_vertices (aBC);
}

/**
 * Sets the values at the outer vertices near the interface (for one type of the elements).
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::sum_up_near_if_outer_values_for_
(
	vector_type & u, ///< the vector where to set
	size_t fct, ///< function index to process
	const DoFDistribution * dd, ///< dof distribution of the grid function to reset
	number time, ///< the physical time
	Grid::VertexAttachmentAccessor<ANumber> & aaBC, ///< attachment for the assembled values
	Grid::VertexAttachmentAccessor<AUInt> & aaNumElem ///< attachment accessor for number of visits
)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator t_elem_iterator;
	
//	get the position accessor (here only available in the LSF)
	position_accessor_type & aaPos = m_spLSF->domain()->position_accessor ();

//	grid level of the dof distribution, and the correct level-set function
	int level = dd->grid_level().level ();
	if (level == GridLevel::TOP)
		level = dd->multi_grid()->top_level ();
	SmartPtr<ls_grid_func_type> spLSF = m_vGLData[level].lsf_on_gl;

//	loop the subsets and the elements in every subset
	std::vector<DoFIndex> multInd (1);
	for (int si = 0; si < dd->num_subsets (); si++)
	{
		t_elem_iterator iterEnd = dd->template end<TElem> (si);
		for (t_elem_iterator iter = dd->template begin<TElem> (si); iter != iterEnd; iter++)
		{
			TElem * pElem = *iter;
		
		//	get the corners of the element
			Vertex * vVertex [TElem::NUM_VERTICES];
			MathVector<dim> vCornerCoords [TElem::NUM_VERTICES];
			for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
				vCornerCoords[co] = aaPos[vVertex[co] = pElem->vertex (co)];
		
		//	check whether we are inside
			if (check_elem_lsf
					(TElem::NUM_VERTICES, pElem, si, level, false, vCornerCoords, time) != 0)
				continue; // this element is not cut, do not consider it
		
		//	get/extrapolate the values of the grid function there
			number vValue [TElem::NUM_VERTICES];
			
		//	get the original values at the corners that are inside
			for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
			{
				if (dd->inner_dof_indices (vVertex[co], fct, multInd) != 1)
					UG_THROW ("LevSetGFsimpleExtrapolation: More than one DoF per vertex for a component. Not the Lagrange element?");
				if (corner_inside (co))
					vValue[co] = DoFRef (u, multInd[0]);
			}
		//	extrapolate values to the corners that are outside and add the contribution
			extrapolate_by_lsf (TElem::NUM_VERTICES, vValue, fct);
			for (size_t co = 0; co < TElem::NUM_VERTICES; co++)
				if (! corner_inside (co))
				{
					aaBC[vVertex[co]] += vValue[co]; // add the value
					(aaNumElem[vVertex[co]])++; // increase the counter
				}
		}
	}
}

/**
 * Sets all matrices at outer vertices to identity.
 */
template <typename TDomain, typename TAlgebra>
void LevSetGFsimpleExtrapolation<TDomain, TAlgebra>::set_outer_matrices
(
	matrix_type & A, ///< the matrix where to set
	const DoFDistribution * dd ///< dof distribution of the grid function to reset
) const
{
	typedef typename DoFDistribution::traits<Vertex>::const_iterator t_vert_iterator;
	
//	If no LSF given, do nothing
	if (m_spLSF.invalid ()) return;

//	Grid level of the dof distribution, and the correct level-set function
	int level = dd->grid_level().level ();
	if (level == GridLevel::TOP)
		level = dd->multi_grid()->top_level ();
	if (level < 0 || (size_t) level >= m_vGLData.size ())
		UG_THROW ("Attempt to assemble on a grid level where the LSF is undefined.");
	SmartPtr<ls_grid_func_type> spLSF = m_vGLData[level].lsf_on_gl;

//	Arrays for the indices in the grid functions:
	std::vector<size_t> vLSFVertInd (1);
	std::vector<size_t> vMatVertInd;
	
//	Loop the vertices
	t_vert_iterator iter = dd->template begin<Vertex> ();
	t_vert_iterator iterEnd = dd->template end<Vertex> ();
	for (; iter != iterEnd; iter++)
	{
		Vertex * pVertex = *iter;
		
	//	Check if the vertex is excluded:
		if (m_excl_ssg.size () != 0 && m_excl_ssg.contains (m_excl_ssg.subset_handler()->get_subset_index (pVertex)))
			continue;
		
	//	Get the multiindex of the LSF and check the value of the LSF
		if (spLSF->inner_algebra_indices (pVertex, vLSFVertInd) != 1)
			UG_THROW ("LevSetGFsimpleExtrapolation: Non-scalar Level-Set Function.");
		if (lsf_inside (BlockRef ((* spLSF) [vLSFVertInd[0]], 0)))
			continue;
		
	//	Get the multiindices of the grid function and set the values:
		size_t n_dofs = dd->inner_algebra_indices (pVertex, vMatVertInd);
		for (size_t dof = 0; dof < n_dofs; dof++)
			SetDirichletRow (A, vMatVertInd[dof]);
	}
}

} // namespace LevelSet
} // end namespace ug

/* End of File */
