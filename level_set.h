/*
 * level_set_util.h
 *
 *  Created on: 01.07.2011
 *      Author: Christian Wehner  christian.wehner@gcsc.uni-frankfurt.de
 */

#ifndef LEVEL_SET_UTIL_H_
#define LEVEL_SET_UTIL_H_

#include "common/common.h"

#include <vector>
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/disc_util/finite_volume_geometry.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include <boost/function.hpp>
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace LevelSet{

template<typename TGridFunction>
class FV1LevelSetDisc
{
	    enum UserDataType {HardcodedData,FunctorData,VectorData,ConstantData};
	///	domain type
		typedef typename TGridFunction::domain_type domain_type;
		
	///	algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

	///	world dimension
		static const int dim = domain_type::dim;

		///for debug	type of scv-size attachment
			//	typedef typename Grid::VertexAttachmentAccessor<Attachment<number> > aaDiv;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	///	type of volume-size attachment
		typedef typename Grid::VertexAttachmentAccessor<Attachment<number> > aaVol;

	///	type of gradient attachment
		typedef typename Grid::VertexAttachmentAccessor<Attachment<MathVector<dim> > > aaGrad;

		// 	Type of multi index vector
		typedef std::vector<MultiIndex<2> > multi_index_vector_type;
		
		// edge iterator
		typedef geometry_traits<EdgeBase>::const_iterator EdgeBaseConstIterator;

		// vertex base iterator
		typedef geometry_traits<VertexBase>::const_iterator VertexBaseConstIterator;

	 public:
    ///	Constructor
      	    ///	Destructor
      	~FV1LevelSetDisc() {};

    ///	set time step
     	FV1LevelSetDisc() :
      		m_dt(0.0), m_time(0), m_gamma(1.0), m_delta(0.0),
      		m_reinit(false), m_analyticalSolution(true),
      		m_analyticalVelocity(false),m_externalVelocity(true),m_analyticalSource(false),m_divFree(false),
      		m_nrOfSteps(1),
      		m_maxCFL(0),m_print(false),m_timestep_nr(0),m_limiter(false),
		    m_vel_x_fct(0),
		    m_vel_y_fct(0),
		    m_vel_z_fct(0),
		    m_source_fct(0),
		    m_solution_fct(0),
		    m_vel_x_vec(0),
		    m_vel_y_vec(0),
		    m_vel_z_vec(0),
		    m_source_vec(0),
		    m_source_constant(0),
		    m_constantv_x(0),
		    m_constantv_y(0),
		    m_constantv_z(0),
		    m_dirichlet_constant(0),
		    m_source_type(HardcodedData),
		    m_velocity_type(HardcodedData),
      	    m_dirichlet_data_type(HardcodedData),
      	    m_exactcurvatureknown(false),
      	    m_interpolate_v_in_ip(true),
      	    m_inside_elements_si(2),
      	  	m_outside_elements_si(3),
      	  	m_onls_elements_si(4),
      	  	m_inside_nodes_si(5),
      	  	m_outside_nodes_si(6),
      	  	m_onls_nodes_si(7)
      	{ set_source(0.0); }

        void set_dt(number deltaT){ UG_LOG("Set dt="<<deltaT<<"\n"); m_dt=deltaT; };

	/// set scale parameters for external velocity and velocity in normal direction
    	void set_vel_scale(number gamma,number delta){ };
    	void set_reinit(size_t n){ m_reinit=1;m_gamma=0;m_delta=1;m_nrOfSteps=n; };
		void set_divfree_bool(bool b){m_divFree=b;};
		void set_gamma(number gamma){m_gamma =gamma;}
		void set_delta(number delta){m_delta =delta;}
		void set_time(double t){m_time = t;}
		number get_time(){return m_time;};
		void set_info(bool b){m_print=b;};
		void set_limiter(bool b){m_limiter=b;};
		void set_timestep_nr(size_t n){m_timestep_nr = n;};
		// set nr of time steps to perform in advect_lsf (default is 1)
		void set_nr_of_steps(size_t n){m_nrOfSteps = n;};
		///	adds a post process to be used when stepping the level set function
//		void add_post_process(IConstraint<dof_distribution_type, algebra_type>& pp){m_vPP.push_back(&pp);}
	    bool compute_error(TGridFunction& numsol);
		bool advect_lsf(TGridFunction& uNew,TGridFunction& u);
	    bool init_function(TGridFunction& u);
	///	adds a post process to be used when stepping the level set function
		void add_post_process(SmartPtr<IConstraint<algebra_type> > pp) {m_vPP.push_back(pp);}

		void set_velocity(SmartPtr<UserData<MathVector<dim>, dim> > user) {m_imVelocity=user;}

		void set_velocity(number vel_x)
		{
			SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
			vel->set_entry(0, vel_x);
			if (dim>1) vel->set_entry(1, vel_x);
			if (dim>2) vel->set_entry(2, vel_x);
			set_velocity(vel);
		}

		void set_velocity(number vel_x, number vel_y)
		{
			SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
			vel->set_entry(0, vel_x);
			vel->set_entry(1, vel_y);
			if (dim>2){
				UG_THROW("ConvectionDiffusion: Setting velocity vector of dimension 2"
									" to a Discretization for world dim " << dim);
			}
			set_velocity(vel);
		}

		void set_velocity(number vel_x, number vel_y, number vel_z)
		{
			SmartPtr<ConstUserVector<dim> > vel(new ConstUserVector<dim>());
			vel->set_entry(0, vel_x);
			vel->set_entry(1, vel_y);
			vel->set_entry(2, vel_z);
			if (dim<3){
				UG_THROW("ConvectionDiffusion: Setting velocity vector of dimension 3"
													" to a Discretization for world dim " << dim);
			}
			set_velocity(vel);
		}

#ifdef UG_FOR_LUA
		void set_velocity(const char* fctName)
		{
			set_velocity(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
		}
#endif

		void set_source(SmartPtr<UserData<number,dim> > user){m_imSource = user;};
		void set_source(number val){set_source(CreateSmartPtr(new ConstUserNumber<dim>(val)));}
#ifdef UG_FOR_LUA
		void set_source(const char* fctName)
		{
			set_source(LuaUserDataFactory<number,dim>::create(fctName));
		}
#endif

		void set_dirichlet_data(SmartPtr<UserData<number,dim> > d){m_imDirichlet = d;};
		void set_dirichlet_data(number val){set_dirichlet_data(CreateSmartPtr(new ConstUserNumber<dim>(val)));};
#ifdef UG_FOR_LUA
		void set_dirichlet_data(const char* fctName)
		{
			set_dirichlet_data(LuaUserDataFactory<number,dim>::create(fctName));
		}
#endif

		bool fill_v_vec(TGridFunction& vel,int component);
		
		bool runtimetest (TGridFunction& u);

		bool compute_normal(TGridFunction& vx,TGridFunction& vy,TGridFunction& u);
		bool compute_dnormal(TGridFunction& dnormal,TGridFunction& vx,TGridFunction& vy,TGridFunction& phi,TGridFunction& u);
		bool compute_ddnormal(TGridFunction& ddnormal,TGridFunction& dnormal,TGridFunction& vx,TGridFunction& vy,TGridFunction& phi,TGridFunction& u);
		
		bool computeElementCurvature2d(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order);
		bool computeElementCurvatureOnGrid2d(TGridFunction& u);

      /// boundary condition subset handling
		bool set_dirichlet_boundary(TGridFunction& uNew,const char* subsets){
			try{
				m_dirichlet_sg = uNew.subset_grp_by_name(subsets);
			}UG_CATCH_THROW("ERROR while parsing Subsets.");

			return true;
       };

	   bool set_outflow_boundary(TGridFunction& uNew,const char* subsets){
			try{
				m_neumann_sg = uNew.subset_grp_by_name(subsets);
			}UG_CATCH_THROW("ERROR while parsing Subsets.");

			return true;
		}

/// subset handling methods:

	   bool init_ls_subsets(TGridFunction& phi);
	   void create_ls_subsets(TGridFunction& phi);
	   bool update_ls_subsets(TGridFunction& phi);

	   void set_elements_active(int sign){
	        if (sign==-1) if (m_inactive_sg.contains(m_inside_elements_si)==true) m_inactive_sg.remove(m_inside_elements_si);
	   		if (sign==0) if (m_inactive_sg.contains(m_onls_elements_si)==true) m_inactive_sg.remove(m_onls_elements_si);
	        if (sign==1) if (m_inactive_sg.contains(m_outside_elements_si)==true)  m_inactive_sg.remove(m_outside_elements_si);
	   	}
	   	void set_elements_active(int signi,int signj){
	       	set_elements_active(signi);
	   	    set_elements_active(signj);
	   	}
	   	void set_elements_active(int signi,int signj,int signk){
	       	set_elements_active(signi,signj);
	   	    set_elements_active(signk);
	   	}
	   	void set_elements_inactive(int sign){
	       	if (sign==-1) m_inactive_sg.add(m_inside_elements_si);
	        if (sign==0) m_inactive_sg.add(m_onls_elements_si);
	        if (sign==1) m_inactive_sg.add(m_outside_elements_si);
	   	}
	   	void set_elements_inactive(int signi,int signj){
	   		set_elements_inactive(signi);
	        set_elements_inactive(signj);
	   	}
	   	void set_elements_inactive(int signi,int signj,int signk){
	       	set_elements_inactive(signi,signj);
	   	    set_elements_inactive(signk);
	   	}
	    void set_nodes_active(int sign){
	   	    if (sign==-1) m_inactive_sg.remove(m_inside_nodes_si);
	      	if (sign==0) m_inactive_sg.remove(m_onls_nodes_si);
	       	if (sign==1) m_inactive_sg.remove(m_outside_nodes_si);
	   	}
	   	void set_nodes_active(int signi,int signj){
	        set_nodes_active(signi);
	        set_nodes_active(signj);
	   	}
	   	void set_nodes_active(int signi,int signj,int signk){
	        set_nodes_active(signi,signj);
	        set_nodes_active(signk);
	   	}
	    void set_nodes_inactive(int sign){
	       	if (sign==-1) m_inactive_sg.add(m_inside_nodes_si);
	       	if (sign==0) m_inactive_sg.add(m_onls_nodes_si);
	        if (sign==1) m_inactive_sg.add(m_outside_nodes_si);
	   	}
	   	void set_nodes_inactive(int signi,int signj){
	       	set_nodes_inactive(signi);
	   	    set_nodes_inactive(signj);
	   	}
	   	void set_nodes_inactive(int signi,int signj,int signk){
	       	set_nodes_inactive(signi,signj);
	       	set_nodes_inactive(signk);
	   	}
	   	bool overwrite(TGridFunction&,TGridFunction&,TGridFunction&,int);
        bool overwrite(TGridFunction&,number,TGridFunction&,int);

	 protected:
	    number analytic_solution(number,MathVector<dim>);
		number analytic_source(number,MathVector<dim>);
	    bool analytic_velocity(MathVector<dim>&,number, MathVector<dim>);

	///	fills the scvVolume attachment for all element types
	    bool calculate_vertex_vol(TGridFunction& u, aaVol& aaVolVolume);

//	    template <typename TElem>
//   	    bool calculate_vertex_grad_vol(grid_type& grid,TGridFunction& u, aaGrad& aaGradient, aaVol& aaVolume );

	    bool calculate_vertex_grad_vol(TGridFunction& u, aaGrad& aaGradient, aaVol& aaVolume );

	    bool calculate_vertex_grad_vol_sign(TGridFunction&, aaGrad& ,aaVol& ,TGridFunction& ,int);

		template <typename TElem>
		bool assemble_element(TElem& elem, DimFV1Geometry<dim>& geo, grid_type& grid,TGridFunction& uNew,const TGridFunction& uOld,aaGrad& aaGradient, aaVol& aaVolume );

//fordebug		template <typename TElem>
//		bool assemble_divergence(TElem& elem,grid_type& grid,TGridFunction& uNew,aaDiv& aaDivergence,aaGrad& aaGradient );

		bool assign_dirichlet(TGridFunction&);
		bool limit_grad(TGridFunction& uOld, aaGrad& aaGradient);

		//bool limit_grad_alpha(TGridFunction& uOld,aaGrad& aaGradient,aaAlpha& aaAlpha);

	private:
	///	vector holding all scheduled post processes
		std::vector<SmartPtr<IConstraint<algebra_type> > > m_vPP;
      	number m_dt;
		number m_time;
    	number m_gamma;
      	number m_delta;
    	bool m_reinit;
    	bool m_analyticalSolution;
	    bool m_analyticalVelocity;
	    bool m_externalVelocity;
		bool m_analyticalSource;
    	bool m_divFree;
    	bool m_exactcurvatureknown;
     	size_t m_nrOfSteps;
		number m_maxCFL;
		bool m_print;
		size_t m_timestep_nr;
		size_t m_limiter;
		SubsetGroup m_neumann_sg;
		SubsetGroup m_dirichlet_sg;
		SubsetGroup m_inactive_sg;
		SmartPtr<UserData<number,dim> > m_vel_x_fct;
		SmartPtr<UserData<number,dim> > m_vel_y_fct;
		SmartPtr<UserData<number,dim> > m_vel_z_fct;
		SmartPtr<UserData<number,dim> > m_source_fct;
		SmartPtr<UserData<number,dim> > m_solution_fct;
		TGridFunction* m_vel_x_vec;
		TGridFunction* m_vel_y_vec;
		TGridFunction* m_vel_z_vec;
		TGridFunction* m_source_vec;
		number m_source_constant;
		number m_constantv_x;
		number m_constantv_y;
		number m_constantv_z;
		number m_dirichlet_constant;
		UserDataType m_source_type;
		UserDataType m_velocity_type;
		UserDataType m_dirichlet_data_type;
		bool m_interpolate_v_in_ip;
		int m_inside_elements_si;
		int m_outside_elements_si;
		int m_onls_elements_si;
		int m_inside_nodes_si;
		int m_outside_nodes_si;
		int m_onls_nodes_si;

		static const size_t maxNumCo = 20;

		///	Data import for the Velocity field
		SmartPtr<UserData<MathVector<dim>, dim> > m_imVelocity;
		///	Data import for the right-hand side
		SmartPtr<UserData<number,dim> > m_imSource;
		///	Data import for the Dirichlet values
		SmartPtr<UserData<number,dim> > m_imDirichlet;
};

} // end namespace LevelSet
} // end namespace ug

// include implementation
#include "level_set_impl.h"

#endif /* LEVEL_SET_UTIL_H_ */
