/*
 * ls_curvature2d.h
 *
 *  Created on: 23.03.2015
 *      Author: Christian Wehner, Dmitriy Logashenko
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET__LS_CURVATURE2D_H__
#define __H__UG__PLUGINS__LEVEL_SET__LS_CURVATURE2D_H__

#include <vector>

#include "common/common.h"

namespace ug{
namespace LevelSet{

/**
 * Analytic functions for computations with the Level-Set-Method
 */
template<typename TGridFunction>
class LevelSetCurvature
{
///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

///	world dimension
	static const int dim = domain_type::dim;

public:

///	class constructor
	LevelSetCurvature()
	:	m_exactcurvatureknown (false)
	{};

///	sets the exact curvature for comparisons
	void exact_curvature(number kappa)
	{
		m_exactcurvatureknown = true;
		m_exactcurv = kappa;
	};

	bool computeElementCurvature2d(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order);
	bool computeElementCurvature2d2(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order,number nodefactor);
	bool computeElementCurvatureOnGrid2d(TGridFunction& u,size_t order,number leastSquaresFactor);
	bool computeElementCurvatureFromSides(TGridFunction& u,size_t order,number leastSquaresFactor);

private:

	bool m_exactcurvatureknown;
	number m_exactcurv;
};

} // end namespace LevelSet
} // end namespace ug

// include implementation
#include "ls_curvature2d_impl.h"

#endif /* __H__UG__PLUGINS__LEVEL_SET__LS_CURVATURE2D_H__ */

/* End of File */
