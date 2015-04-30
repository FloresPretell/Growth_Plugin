/*
 * ls_analytic_impl.h
 *
 *  Created on: 23.03.2015
 *      Author: Christian Wehner, Dmitriy Logashenko
 */

namespace ug{
namespace LevelSet{

template<typename TGridFunction>
bool LevelSetAnalytic<TGridFunction>::analytic_velocity(MathVector<dim> & v,number t,MathVector<dim> x)
{
	switch(dim)
	{
	    case 1: 
		   v[0] = 1;
		   return 0;
		case 2: 
		   v[0] = -x[1];
		   v[1] = x[0];
		  // v[0] = 1;
		  // v[1] = 0;
//		   v[0] = 1;
	//	   v[1] = 1;
		   return 0;
		case 3: 
		   v[0] = -x[1];
		   v[1] =  x[0];
		   v[2] =  0;
		   return 0;
	};		
	return true;
};

template<typename TGridFunction>
number LevelSetAnalytic<TGridFunction>::analytic_source(number t,MathVector<dim> x)
{
	switch(dim)
	{
	    case 1: 
		   return 1;
		case 2: 
		   return 0;
		case 3: 
		   return 1;
	};		
	return 1;
};

template<typename TGridFunction>
number LevelSetAnalytic<TGridFunction>::analytic_solution(number t,MathVector<dim> x)
{
	MathVector<dim> xnew;
    switch(dim)
	{
        case 1:  
			//number v=1;
			//xnew=x-v*t;
			//return sin(xnew);
        	return 0;
        case 2:
        	// return 0;
        	//return x[1];//-t;
          //  number z;
      //      z=min(min(std::abs(x[0]-1),std::abs(x[0]+1)),min(std::abs(x[1]-1),std::abs(x[1]+1)));
            //z=min(min(std::abs(x[0]),std::abs(1-x[0])),min(std::abs(x[1]),std::abs(1-x[1])));
        //	return z;
		    xnew[0] = x[0]*cos(t)+x[1]*sin(t);
		 	xnew[1] = -x[0]*sin(t)+x[1]*cos(t);
		 	return sqrt((xnew[0]-0.5)*(xnew[0]-0.5) + xnew[1]*xnew[1]) - 0.5;
            xnew[0] = x[0] - t;
            xnew[1] = x[1] - t;
			return x[0]-t;
			return sqrt(xnew[0]*xnew[0] + xnew[1]*xnew[1]) - 0.5;
           // return 2*xnew[0]-xnew[1];
        case 3:  
			return 0;
	}

    return -1;
};

template<typename TGridFunction>
bool LevelSetAnalytic<TGridFunction>::fill_v_vec(TGridFunction& vel,int component)
{
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexConstIterator;

	if (component>dim){
		UG_THROW("fill_v_vec: component > dim.");
	}

	position_accessor_type aaPos = vel.domain()->position_accessor();
	for (int si=0;si<vel.num_subsets();++si)
	{
		for(VertexConstIterator iter = vel.template begin<Vertex>(si);
									   iter != vel.template end<Vertex>(si); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
			MathVector<dim> coord;
			MathVector<dim> vnode;
			coord = aaPos[vrt];

		//	get vector holding all indices on the vertex
			std::vector<DoFIndex> ind;

			vel.inner_dof_indices(vrt, 0, ind);
			analytic_velocity(vnode,m_time,coord);
			DoFRef(vel, ind[0]) = vnode[component];
	     }
	};
	return true;
};

// initialize level set function with analytical solution
template<typename TGridFunction>
bool LevelSetAnalytic<TGridFunction>::init_function(TGridFunction& u)
{
//	get domain of grid function
	domain_type& domain = *u.domain().get();

//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	typedef typename domain_type::position_accessor_type position_accessor_type;

	//	read indices on vertex
	position_accessor_type aaPos = domain.position_accessor();

	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexConstIterator;
	for (int si=0;si<domain.subset_handler()->num_subsets();++si){
		for(VertexConstIterator iter = u.template begin<Vertex>(si);
									   iter != u.template end<Vertex>(si); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
		//	get vector holding all indices on the vertex
			std::vector<DoFIndex> ind;
			u.inner_dof_indices(vrt, 0, ind);
			DoFRef(u, ind[0]) = analytic_solution(m_time,coord);
	     }
	};
	return true;
};

} // end namespace LevelSet
} // end namespace ug

/* End of File */
