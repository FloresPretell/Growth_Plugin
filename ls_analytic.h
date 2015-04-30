/*
 * ls_analytic.h
 *
 *  Created on: 23.03.2015
 *      Author: Christian Wehner, Dmitriy Logashenko
 */
#ifndef __H__UG__PLUGINS__LEVEL_SET__LS_ANALYTIC_H__
#define __H__UG__PLUGINS__LEVEL_SET__LS_ANALYTIC_H__

#include "common/common.h"

namespace ug{
namespace LevelSet{

/**
 * Analytic functions for computations with the Level-Set-Method
 */
template<typename TGridFunction>
class LevelSetAnalytic
{
///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

///	world dimension
	static const int dim = domain_type::dim;

public:

///	class constructor
	LevelSetAnalytic()
	:	m_time (0)
	{};

	void set_time(number t) {m_time = t;}
	number get_time() {return m_time;};
	
	number analytic_solution(number,MathVector<dim>);
	number analytic_source(number,MathVector<dim>);
	bool analytic_velocity(MathVector<dim>&,number, MathVector<dim>);

	bool fill_v_vec(TGridFunction& vel,int component);
	bool init_function(TGridFunction& u);
	
private:

	number m_time;
};

} // end namespace LevelSet
} // end namespace ug

// include implementation
#include "ls_analytic_impl.h"

#endif /* __H__UG__PLUGINS__LEVEL_SET__LS_ANALYTIC_H__ */

/* End of File */
