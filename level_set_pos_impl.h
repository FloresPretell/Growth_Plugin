/*
 * Tools for determining the position of the level set: Implementation
 *
 * Created on Jul. 10, 2015 by D. Logashenko
 */

namespace ug{
namespace LevelSet{

/**
 * get the sign of the LSF in an element:
 * 0 if there both the signes,
 * -1 if negative at all the corners
 * 1 if positive at all the corners
 */
template<typename TGridFunction>
int LSPositionZ<TGridFunction>::lsf_sign
(
	size_t noc, ///< number of corners
	number lsf [] ///< corner values of the LSF
)
{
	int pos = 0, neg = 0;
	for (size_t co = 0; co < noc; co++)
	{
		if (lsf [co] >= - lsf_threshold () && lsf [co] <= lsf_threshold ())
			return 0; // we consider such elements as intersected; this is important for computation of the gradients
		if (lsf [co] < 0) neg = 1;
		else pos = 1;
	}
	return pos - neg;
}

///	Find the cut elements
template <typename TGridFunction>
void LSPositionZ<TGridFunction>::fill_cut ()
{
///	grid element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator t_elem_iter;

//	reset the array
	m_vpCut.clear ();
	
//	if no LSF then nothing to do
	if (! m_spLSF.valid ()) return;
	
	LocalIndices locInd;
	LocalVector locLSF;
	
//	loop the vertices
	t_elem_iter iterEnd = m_spLSF->template end<elem_type> ();
	for (t_elem_iter iter = m_spLSF->template begin<elem_type> (); iter != iterEnd; ++iter)
	{
		elem_type * elem = *iter;
		size_t noc = elem->num_vertices ();
		number lsf [maxNumCorners];
		
	//	local values of the LSF
		m_spLSF->indices (elem, locInd);
		locLSF.resize (locInd);
		GetLocalVector (locLSF, *m_spLSF);
		for (size_t i = 0; i < noc; i++) lsf[i] = locLSF (0, i);
		
	//	check the sign of the lsf
		if (lsf_sign (noc, lsf) == 0)
			m_vpCut.push_back (elem);
	}
}

///	Find the elements intersected by the straight line
template <typename TGridFunction>
bool LSPositionZ<TGridFunction>::get_intersected
(
	MathVector<dim-1> xy ///< the first coordinates of the points of the line
)
{
	MathVector<dim> origin, dir;
	
//	prepare the full-dimensional vectors
	VecSet (dir, 0); dir[dim-1] = 1;
	for (size_t i = 0; i < dim-1; i++) origin[i] = xy[i];
	origin[dim-1] = 0;
	
//	get the intersections
	m_intersections.clear ();
	if (! RayElementIntersections (m_intersectionRecords, m_tree, origin, dir))
		return false;
	for (size_t i = 0; i < m_intersectionRecords.size (); i++)
	{
		intersect_rec_type & r = m_intersectionRecords[i];
		m_intersections.push_back (elem_intersect_data (r.elem,
			PointOnRay(origin, dir, r.smin), PointOnRay(origin, dir, r.smax)));
	}
	return true;
}

///	Compute the coordinates of the intersection with the level set
template <typename TGridFunction>
void LSPositionZ<TGridFunction>::get_z
(
	std::vector<number> & z, ///< the z-coordinates of the intersections
	number small_z ///< a threshold for the coordinates
)
{
	LocalIndices locInd;
	LocalVector locLSF;
	std::vector<number> vShape;
		
	for (size_t i = 0; i < m_intersections.size (); i++)
	{
		elem_intersect_data & e_data = m_intersections[i];
		elem_type * elem = e_data.elem;
		
	//	get corners of element
		std::vector<MathVector<dim> > vCornerCoords;
		CollectCornerCoordinates (vCornerCoords, *elem, *m_spLSF->domain ());

	//	reference object id, the reference mapping and the trial space
		const ReferenceObjectID roid = elem->reference_object_id ();
		DimReferenceMapping<dim, dim> & map
			= ReferenceMappingProvider::get<dim, dim> (roid, vCornerCoords);
		const LocalShapeFunctionSet<dim>& rTrialSpace =
				LocalFiniteElementProvider::get<dim> (roid, m_spLSF->local_finite_element_id (0));

	//	get local position of DoF
		MathVector<dim> locPos [2];
		VecSet (locPos[0], 0.5); VecSet (locPos[1], 0.5);
		map.global_to_local(locPos, e_data.pnt, 2);

	//	local values of the LSF
		m_spLSF->indices (elem, locInd);
		locLSF.resize (locInd);
		GetLocalVector (locLSF, *m_spLSF);

	//	process the points
		number lsf_value [2];
		for (size_t j = 0; j < 2; j++)
		{
		//	evaluate at shapes at the point
			rTrialSpace.shapes (vShape, locPos[j]);

		// 	compute value of the LSF at the point
			lsf_value[j] = 0.0;
			for(size_t sh = 0; sh < vShape.size(); sh++)
				lsf_value[j] += locLSF (0, sh) * vShape[sh];
		}
		
	//	check if the element contains the intersection of the level set with the line
		if (lsf_value[0] * lsf_value[1] > 0) // stricktly ">" not to loose the intersections at faces
			continue;
		
	//	find the intersection (assume the linearity of the LSF along the line)
		number elem_z = (lsf_value[0] * e_data.pnt[0][dim-1] - lsf_value[1] * e_data.pnt[1][dim-1])
						/ (lsf_value[0] - lsf_value[1]);
		
	//	check if we have already this value
		size_t k;
		for (k = 0; k < z.size (); i++)
			if (std::fabs (z[k] - elem_z) < small_z) break;
		if (k >= z.size ()) // i.e. if none found
			z.push_back (elem_z);
	}
}

} // end namespace LevelSet
} // end namespace ug

/* End of File */
