/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Authors: Christian Wehner, Dmitriy Logashenko
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
