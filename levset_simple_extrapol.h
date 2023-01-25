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
 * Extrapolation over the interface based on the simple linear extrapolation along the ray.
 */
#ifndef __H__UG__PLUGINS__LEVSET_SIMPLE_EXTRAPOL__
#define __H__UG__PLUGINS__LEVSET_SIMPLE_EXTRAPOL__

#include <vector>

// ug4 headers
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/dom_disc_embb.h"

namespace ug {
namespace LevelSet {

/// Level-set extrapolation for the ghost-fluid method based on piecewise linear LS functions
/**
 * This class implements the extrapolation of the vertex-centered solution
 * over the interface given by the zero isoline of a piecewise linear vertex-centered
 * grid function.
 *
 * \tparam TDomain	domain type
 * \tparam TAlgebra	algebra type for the functions to extrapolate
 */
template <typename TDomain, typename TAlgebra>
class LevSetGFsimpleExtrapolation
{
	typedef LevSetGFsimpleExtrapolation<TDomain, TAlgebra> this_type;
	
public:

///	the threshold for the level-set function
	static constexpr number small_lsf_value = 1e-8;
	
///	returns true if the given value should correspond to the 'inside' subdomain
	/**
	 * Note that this function is used only to determine whether a vertex is
	 * inside the actual domain. It introduces a threshold that should exclude
	 * machine zero in denominators in the extrapolation and the elimination in
	 * the Ghost-Fluid method. In the extrapolation and the elimination themselves,
	 * no threshold is used, and the interface corresponds to the zero isoline
	 * of the LSF.
	 */
	static bool lsf_inside
	(
		number lsf_val ///< value of the level-set function
	)
	{
		return lsf_val < - small_lsf_value;
	}
	
///	domain type
	typedef TDomain domain_type;
	
///	grid type for the domain
	typedef typename TDomain::grid_type grid_type;
	
///	algebra type for the functions to extrapolate
	typedef TAlgebra algebra_type;
	
///	algebra type for the LSF
	typedef CPUAlgebra lsf_algebra_type;
	
///	vector type (for the functions to extrapolate)
	typedef typename algebra_type::vector_type vector_type;
	
///	matrix type
	typedef typename algebra_type::matrix_type matrix_type;
	
///	grid function type for the LSF
	typedef GridFunction<domain_type, lsf_algebra_type> ls_grid_func_type;
	
///	type of the projection of the LSF to the coarser grid levels
	typedef StdInjection<domain_type, lsf_algebra_type> projection_type;
	
///	type of approximation space
	typedef ApproximationSpace<domain_type> approx_space_type;
	
///	type of the position attachment accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	
///	dimensionality (the World dimension)
	static const int dim = domain_type::dim;
	
private:

///	types of boundary conditions at the moving boundary
	enum t_fs_bc
	{
		DIRICHLET_BC = 0,
		DIRICHLET_PLAIN_BC,
		NEUMANN_0_BC,
		MAX_FS_BC
	};
	
///	max. number of corners of an grid element
	static const size_t max_num_corners = grid_dim_traits<dim>::MaxNumVerticesOfElem;
	
public:

/// class constructor
	LevSetGFsimpleExtrapolation ()
	:	m_spLSF (SPNULL)
	{}
	
///	set the level-set function and check it
	void set_LSF
	(
		SmartPtr<ls_grid_func_type> spLSF ///< the function to set
	)
	{
		if (spLSF.valid () && (spLSF->local_finite_element_id (0) != LFEID (LFEID::LAGRANGE, dim, 1)
								|| spLSF->dd()->num_fct () != 1))
			UG_THROW ("LevSetGFExtrapolation: The Level-Set Function must be scalar piecewise linear.");
		m_spLSF = spLSF;
		prepare_grid_levels ();
	}
	
///	projects the values of the LSF to the coarser levels
	void project_LSF ();
	
///	checks the dof distribution of the solution
	void check_dd
	(
		ConstSmartPtr<DoFDistribution> dd ///< the dof distribution to check
	) const
	{
		if (m_spLSF.invalid ()) return; // there is no LSF, so no matter!
		
		if (m_vICData.size () != dd->num_fct ())
			UG_THROW ("LevSetGFExtrapolation: Mismatch of number of the functions!");
		for (size_t fct = 0; fct < dd->num_fct (); fct++)
			if (dd->local_finite_element_id (fct) != LFEID (LFEID::LAGRANGE, dim, 1))
				UG_THROW ("LevSetGFExtrapolation: The Level-Set Function works currently only with piecewise linear functions.");
	}
	
///	prepares the boundary conditions at the interface: sets all them to Dirichlet-0
	void prepare_interface_bc
	(
		SmartPtr<approx_space_type> spApproxSpace ///< the approximation space of the domain discretization
	)
	{
		m_vICData.resize (spApproxSpace->num_fct ());
	}
	
///	adds a Dirichlet BC with a given value
	void set_Dirichlet_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		number value ///< the Dirichlet value
	)
	{
		prepare_interface_bc (spApproxSpace);
		size_t fct = spApproxSpace->fct_id_by_name (fct_name);
		m_vICData[fct].bc_kind = DIRICHLET_BC;
		m_vICData[fct].func = SPNULL; // we do not use the function here
		m_vICData[fct].value = value;
	}
	
///	adds a Dirichlet BC with a given value
	void set_Dirichlet_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		SmartPtr<CplUserData<number, dim> > func ///< the Dirichlet function
	)
	{
		prepare_interface_bc (spApproxSpace);
		size_t fct = spApproxSpace->fct_id_by_name (fct_name);
		m_vICData[fct].bc_kind = DIRICHLET_BC;
		m_vICData[fct].func = func;
		m_vICData[fct].value = 0; // a dummy value
	}
	
///	adds a Dirichlet BC with a given value
	void set_plain_Dirichlet_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		number value ///< the Dirichlet value
	)
	{
		prepare_interface_bc (spApproxSpace);
		size_t fct = spApproxSpace->fct_id_by_name (fct_name);
		m_vICData[fct].bc_kind = DIRICHLET_PLAIN_BC;
		m_vICData[fct].func = SPNULL; // we do not use the function here
		m_vICData[fct].value = value;
	}
	
///	adds a Dirichlet BC with a given value
	void set_plain_Dirichlet_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		SmartPtr<CplUserData<number, dim> > func ///< the Dirichlet function
	)
	{
		prepare_interface_bc (spApproxSpace);
		size_t fct = spApproxSpace->fct_id_by_name (fct_name);
		m_vICData[fct].bc_kind = DIRICHLET_PLAIN_BC;
		m_vICData[fct].func = func;
		m_vICData[fct].value = 0; // a dummy value
	}
	
///	adds a Neumann-0
	void set_Neumann0_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name ///< function to impose the condition for
	)
	{
		prepare_interface_bc (spApproxSpace);
		size_t fct = spApproxSpace->fct_id_by_name (fct_name);
		m_vICData[fct].bc_kind = NEUMANN_0_BC;
		m_vICData[fct].func = SPNULL; // a dummy value
		m_vICData[fct].value = 0; // a dummy value
	}
	
///	excludes a (boundary) subsets from the extrapolation; note that this resets the old excluded subsets
	void exclude_subsets
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* subset_names ///< names of the subsets to exclude
	)
	{
		std::string ss_names (subset_names);
		std::vector<std::string> v_ss_names;
		TokenizeTrimString (ss_names, v_ss_names);
		m_excl_ssg.set_subset_handler (spApproxSpace->subset_handler ());
		m_excl_ssg.add (v_ss_names);
	}
	
private:

///	gets the BC value
	number boundary_value
	(
		size_t fct ///< the function
	) const;

///	checks if the corner is excluded
	bool corner_excluded
	(
		size_t co ///< the corner
	) const
	{
		return m_co_excluded [co];
	}

//----	Extrapolation	----//
public:
	
///	checks whether the element is intersected by the interface or not, and prepares the data
	int check_elem_lsf
	(
		size_t n_co, ///< number of the corners of the element
		GridObject * pElem, ///< the element to process
		int si, ///< subset of the element
		int g_level, ///< grid level of the element
		bool use_hanging, ///< if there can be hanging nodes
		const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
		number time ///< the phisical time
	);
	
///	returns true if the corner is "inside" (use after check_elem_lsf)
	bool corner_inside
	(
		size_t co ///< the corner
	) const
	{return lsf_inside (m_locLSF.value (0, co));}
	
///	returns the effective value of the LSF at a corner (use after check_elem_lsf)
	number lsf_at
	(
		size_t co ///< the corner
	) const
	{return m_locLSF.value (0, co);}
	
///	extrapolates all the components of the solution to the vertices behind the interface (w.r.t. a base corner)
	template <typename TElem>
	void extrapolate_sol_by_lsf
	(
		LocalVector& locU, ///< the solution to extrapolate
		size_t base_co ///< the base corner of the element
	) const;
	
/// extrapolates a component of the solution to the vertices behind the interface (w.r.t. a base corner)
	void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		size_t base_co, ///< the base corner of the element
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const;
	
///	extrapolates a component of the solution to the vertices behind the interface (by averaging)
	void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const;

///	sets the values of the vector at the non-base vertices to 0
	template <typename TElem>
	void clear_outer_vectors
	(
		LocalVector& locD, ///< the local defect/rhs to postprocess
		size_t base_co ///< the base corner of the element
	);
	
///	eliminates the matrix connections to the vertices behind the interface
	template <typename TElem>
	void eliminate_extrapolated
	(
		LocalMatrix& locM, ///< the local matrix to postprocess
		size_t base_co ///< the base corner of the element
	);

///	eliminates the matrix connections to the vertices behind the interface and the corresponding rhs
	template <typename TElem>
	void eliminate_extrapolated
	(
		LocalMatrix& locM, ///< the local matrix to postprocess
		LocalVector& locB, ///< the local right-hand side to postprocess
		size_t base_co ///< the base corner of the element
	);

///	eliminates the matrix connections to the vertices behind the interface and the corresponding rhs
	template <typename TElem>
	void eliminate_extrapolated
	(
		LocalMatrix& locM, ///< the local matrix to postprocess
		LocalVector& locB ///< the local right-hand side to postprocess
	);

//----	Resetting data in the outer parts	----//
public:

	/// sets the values at the outer vertices to 0
	void clear_outer_values
	(
		vector_type & d, ///< the vector where to set
		const DoFDistribution * dd ///< dof distribution of the grid function to reset
	) const;
	
	/// sets the values at the outer vertices to given values
	void set_outer_values
	(
		vector_type & u, ///< the vector where to set
		const DoFDistribution * dd, ///< dof distribution of the grid function to reset
		number time ///< the physical time
	);
	
	/// sets the matrices at outer vertices to identity
	void set_outer_matrices
	(
		matrix_type & A, ///< the vector where to set
		const DoFDistribution * dd ///< dof distribution of the grid function to reset
	) const;
	
private:

	///	sets the values at the outer vertices near the interface (for one type of the elements)
	template <typename TElem>
	void sum_up_near_if_outer_values_for_
	(
		vector_type & u, ///< the vector where to set
		size_t fct, ///< the function index to process
		const DoFDistribution * dd, ///< dof distribution of the grid function to reset
		number time, ///< the physical time
		Grid::VertexAttachmentAccessor<ANumber> & aaBC, ///< attachment accessor for the bc values
		Grid::VertexAttachmentAccessor<AUInt> & aaNumElem ///< attachment accessor for number of visits
	);
	
	/// helper class for the loop over all the element types
	struct SumUpNearIfOuterValues
	{
		this_type * m_pThis;
		vector_type & m_u;
		size_t m_fct;
		const DoFDistribution * m_dd;
		number m_time;
		Grid::VertexAttachmentAccessor<ANumber> & m_aaBC;
		Grid::VertexAttachmentAccessor<AUInt> & m_aaNumElem;
		
		SumUpNearIfOuterValues
		(
			this_type * pThis, ///< 'this' pointer
			vector_type & u, ///< the vector where to set
			size_t fct, ///< the function index to process
			const DoFDistribution * dd, ///< dof distribution of the grid function to reset
			number time, ///< the physical time
			Grid::VertexAttachmentAccessor<ANumber> & aaBC, ///< attachment accessor for the bc values
			Grid::VertexAttachmentAccessor<AUInt> & aaNumElem ///< attachment accessor for number of visits
		) : m_pThis (pThis), m_u (u), m_fct (fct), m_dd (dd), m_time (time), m_aaBC (aaBC), m_aaNumElem (aaNumElem) {};
		
		template <typename TElem> void operator () (TElem &)
		{
			m_pThis->template sum_up_near_if_outer_values_for_<TElem> (m_u, m_fct, m_dd, m_time, m_aaBC, m_aaNumElem);
		}
	};

	///	projects the level-set functions to the coarser grid levels
	void prepare_grid_levels ();
	
private:

///	Class for the level-set function on the coarser levels
	struct GLData
	{
		SmartPtr<ls_grid_func_type> lsf_on_gl; ///< level-set function on the level
		SmartPtr<projection_type> inject; ///< the canonical restriction for the grid level
	};
	
	SmartPtr<ls_grid_func_type> m_spLSF; ///< original Level-Set function
	std::vector<GLData> m_vGLData; ///< data on the grid levels
	
///	Class for the specification of the boundary conditions at the interface (for one component/function)
	struct ICData
	{
		t_fs_bc bc_kind; ///< type of the interface conditions on the free surface
		SmartPtr<CplUserData<number, dim> > func; ///< variable value (for the BC)
		number value; ///< value (for the Dirichlet BC, if func is NULL)
		
		/// Default constructor
		ICData () : bc_kind (DIRICHLET_BC), func (SPNULL), value (0) {} // dummy values
	};
	
	std::vector<ICData> m_vICData; ///< the boundary conditions at the interface
	
	SubsetGroup m_excl_ssg; ///< subsets to exclude
	
//	temporaty data (valid for the initialization with a particular element)

	LocalIndices m_indLSF; ///< indices of the LSF dof's
	LocalVector m_locLSF; ///< corner values of the LFS (if initialized)

	GridObject * m_pElem; ///< the current grid element
	int m_si; ///< subset index of the element
	MathVector<dim> m_elemCenter; ///< center of the element (only initialized for cut elements)
	const MathVector<dim> * m_vCornerCoords; ///< array of coordinates of the corners
	bool m_co_excluded [max_num_corners]; ///< if corners are excluded by subsets (only for cut elements)
	
	number m_time; ///< the physical time from the current initialization
	
}; // class LevSetGFExtrapolation

/// Global assembler for the simple extrapolation method
/**
 * Global assembler for the simple extrapolation method based on the LevSetGFsimpleExtrapolation
 * class.
 *
 * \tparam TDomain		domain type
 * \tparam TAlgebra		algebra type
 */
template <typename TDomain, typename TAlgebra>
class LSGFsimpleDomainDiscretization
:	public LSGFDomainDiscretization<TDomain, TAlgebra, LevSetGFsimpleExtrapolation<TDomain, TAlgebra> >
{
///	Type of the base class
	typedef LSGFDomainDiscretization<TDomain, TAlgebra, LevSetGFsimpleExtrapolation<TDomain, TAlgebra> > base_type;
	
///	Type of approximation space
	typedef ApproximationSpace<TDomain>	approx_space_type;
	
public:

///	default Constructor
	LSGFsimpleDomainDiscretization (SmartPtr<approx_space_type> pApproxSpace)
	: 	base_type (pApproxSpace)
	{}
};

} // namespace LevelSet
} // end namespace ug

#include "levset_simple_extrapol_impl.h"

#endif // __H__UG__PLUGINS__LEVSET_SIMPLE_EXTRAPOL__

/* End of File */
