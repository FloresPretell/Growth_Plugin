/*
 * level_set_user_data.h
 *
 *  Created on: 15.11.2012
 *      Author: Christian Wehner
 */

#ifndef __LEVEL_SET__LEVEL_SET_USER_DATA__
#define __LEVEL_SET__LEVEL_SET_USER_DATA__

#include "common/common.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/lib_disc.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"


namespace ug{
namespace LevelSet{

template <typename TData, int dim, typename TImpl>
class LevelSetUserDataBase
: 	public StdUserData<LevelSetUserDataBase<TData,dim,TImpl>, TData, dim>
{
	public:
	virtual void operator() (TData& value,
	                         const MathVector<dim>& globIP,
	                         number time, int si) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void operator() (TData vValue[],
	                         const MathVector<dim> vGlobIP[],
	                         number time, int si, const size_t nip) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	template <int refDim>
	void evaluate(TData vValue[],
	                const MathVector<dim> vGlobIP[],
	                number time, int si,
	                GridObject* elem,
	                const MathVector<dim> vCornerCoords[],
	                const MathVector<refDim> vLocIP[],
	                const size_t nip,
	                LocalVector* u,
	                const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		getImpl().template evaluate<refDim>(vValue,vGlobIP,time,si,elem,
		                                    vCornerCoords,vLocIP,nip,u,vJT);
	}

	virtual void compute(LocalVector* u, GridObject* elem,
	                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const number t = this->time();
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			getImpl().template evaluate<dim>(this->values(s), this->ips(s), t, si,
			                                 elem, vCornerCoords, this->template local_ips<dim>(s),
			                                 this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
	                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			getImpl().template evaluate<dim>(this->values(s), this->ips(s),  this->time(s), si,
			                                 elem, vCornerCoords, this->template local_ips<dim>(s),
			                                 this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///	returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///	returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}

	protected:
	///	access to implementation
	TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
	const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};



template <typename TGridFunction>
class LevelSetUserData
: public LevelSetUserDataBase<number, TGridFunction::dim,
  LevelSetUserData<TGridFunction> >, virtual public INewtonUpdate
  {
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
	static const int dim = domain_type::dim;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	enum eval_type {sharp,cr_ip_average};

	  public:
	// set evaluation type, implemented so far:
	// 0 sharp (compute lsf in ip and give out inside value if inside or outside value if outside)
	// 1 cr_ip_average (compute averaged value in ip, if in CR-FV-Geometry from-value in ip and to-value in ip have different signs,
	//                  harmonic average is computed using intersection position from from and to node which is computed from level set function,
	//                  the ips given in evaluate must be the CR-FV ips)
	void set_eval_type(int type){
		if (type==0) m_eval_type=sharp;
		if (type==1) m_eval_type=cr_ip_average;
	}

	eval_type m_eval_type;

	  private:
	// level set grid function
	SmartPtr<TGridFunction> m_phi;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

	/*	//	component of function
		size_t m_fct;

	//	local finite element id
		LFEID m_lfeID;               */

	  public:
	void set_inside_data(SmartPtr<CplUserData<number, dim> > user){
		m_imInsideData = user;
	}
	void set_inside_data(number val){
		set_inside_data(make_sp(new ConstUserNumber<dim>(val)));
	}
#ifdef UG_FOR_LUA
	void set_inside_data(const char* fctName){
		set_inside_data(LuaUserDataFactory<number, dim>::create(fctName));
	}
#endif
	void set_outside_data(SmartPtr<CplUserData<number, dim> > user){
		m_imOutsideData = user;
	}
	void set_outside_data(number val){
		set_outside_data(make_sp(new ConstUserNumber<dim>(val)));
	}
#ifdef UG_FOR_LUA
	void set_outside_data(const char* fctName){
		set_outside_data(LuaUserDataFactory<number, dim>::create(fctName));
	}
#endif

	  private:
	///	Data import for inside and outside data
	SmartPtr<CplUserData<number,dim> > m_imInsideData;
	SmartPtr<CplUserData<number,dim> > m_imOutsideData;

	  public:
	/// constructor
	LevelSetUserData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		m_phi = spGridFct;
		domain_type& domain = *m_phi->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		m_eval_type = cr_ip_average;
	}

	virtual ~LevelSetUserData(){};

	template <int refDim>
	inline void evaluate(number vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GridObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		const size_t numVertices = element->num_vertices();
		std::vector<Vertex*> childVertex(numVertices);
		std::vector<number> phi(numVertices);

	//	for (size_t i=0;i<nip;i++){
	//		UG_LOG("co(" << i << ",:)=" << vGlobIP[i] << "\n");
	//	}

		// find child vertices by injection
		for(size_t i = 0; i < numVertices; ++i){
			childVertex[i] = element->vertex(i);
		};

		// find out level, then get child vertices on finest level
		size_t numChildren = m_grid->template num_children<elem_type>(element);
		if (numChildren!=0){
			size_t lowerLevel=0;
			elem_type* childElem = m_grid->template get_child<elem_type>(element,0);
			numChildren = m_grid->template num_children<elem_type>(childElem);
			lowerLevel++;
			while (numChildren>0){
				childElem = m_grid->template get_child<elem_type>(childElem,0);
				numChildren = m_grid->template num_children<elem_type>(childElem);
				lowerLevel++;
			};
			for (size_t i=0;i<nip;i++){
				for (size_t j=0;j<lowerLevel;j++){
					childVertex[i]=m_grid->template get_child<Vertex>(childVertex[i],0);
				}
			}
		}
		//    create Multiindex
		std::vector<DoFIndex> ind;
		for (size_t i=0;i<numVertices;i++){
			m_phi->dof_indices(childVertex[i], 0, ind);
			phi[i]=DoFRef(*m_phi, ind[0]);
		};
		bool onls=false;
		bool inside=false;
		for (size_t i=0;i<numVertices;i++){
			if (phi[i]==0){
				continue;
			};
			if (phi[i]<0) inside=true;
			for (size_t j=i+1;j<numVertices;j++){
				if (phi[i]*phi[j]<0){
					onls = true;
					break;
				};
			}
		}
		if (onls==false){
			if (inside==false){
				(*m_imInsideData)(vValue,
						vGlobIP,
						time, si,
						elem,
						vCornerCoords,
						vLocIP,
						nip,
						u,
						vJT);
			} else {
				(*m_imOutsideData)(vValue,
						vGlobIP,
						time, si,
						elem,
						vCornerCoords,
						vLocIP,
						nip,
						u,
						vJT);
			};
			return;
		};
		number vValueInside[max_number_of_ips];
		number vValueOutside[max_number_of_ips];
		(*m_imInsideData)(vValueInside,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		(*m_imOutsideData)(vValueOutside,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		if (m_eval_type == sharp){
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	reference object id
				ReferenceObjectID roid = elem->reference_object_id();

				//	memory for shapes
				std::vector<number> vShape;

				// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
						LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

				//  evaluate shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

				//    get multiindices of element
				std::vector<DoFIndex> ind;
				m_phi->dof_indices(elem, 0, ind);

				//     compute lsf at integration point
				number phiValue = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					phiValue += phi[sh] * vShape[sh];
				}
				if (phiValue<0) vValue[ip] = vValueInside[ip];
				else  vValue[ip] = vValueOutside[ip];
				//		                    UG_LOG("vValue(" << ip << ")=" << vValueOutside[ip] << "\n");
			};
		}

		if (m_eval_type == cr_ip_average){
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	memory for shapes
			std::vector<number> vShape;

			//    get domain of grid function
			const domain_type& domain = *m_phi->domain().get();

			//    get position accessor
			typedef typename domain_type::position_accessor_type position_accessor_type;
			const position_accessor_type& posAcc = domain.position_accessor();

			position_accessor_type aaPos = m_phi->domain()->position_accessor();

			//    coord and vertex array
			MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
			Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
			DimCRFVGeometry<dim> geo;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = element->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			};
			//    evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			std::vector<number> phiSideValue(geo.num_sh());
			// compute interpolated lsf values in cr dofs
			for (size_t i=0;i<geo.num_scv();i++){
				// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
				const LocalShapeFunctionSet<dim>& rTrialSpace =
						LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

				//  evaluate shapes at ip
				rTrialSpace.shapes(vShape, geo.scv(i).local_ip());

				//    get multiindices of element
				std::vector<DoFIndex> ind;
				m_phi->dof_indices(elem, 0, ind);

				//     compute lsf at integration point
				phiSideValue[i] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					phiSideValue[i] += phi[sh] * vShape[sh];
				}
			}
			for (size_t ip=0;ip<nip;ip++){
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				number phiFrom = phiSideValue[scvf.from()];
				number phiTo = phiSideValue[scvf.to()];
				//		                    UG_LOG("ip=" << ip << " inside=" << vValueInside[ip] << " outside=" << vValueOutside[ip] << "\n");
				if (phiFrom==0){
					if (phiTo<=0){
						vValue[ip] = vValueInside[ip];
					} else {
						vValue[ip] = vValueOutside[ip];
					}
				} else {
					number theta;
					if (phiFrom<0){
						if (phiTo<=0){
							vValue[ip] = vValueInside[ip];
						} else {
							theta = phiFrom/(phiFrom-phiTo);
							vValue[ip] = vValueInside[ip]*vValueOutside[ip]/(theta*vValueInside[ip]+(1-theta)*vValueOutside[ip]);
							//		                                UG_LOG("theta=" << theta << " vValue=" << vValue[ip] << "\n");
						}
					} else {
						if (phiTo>=0){
							vValue[ip] = vValueOutside[ip];
						} else {
							theta = phiFrom/(phiFrom-phiTo);
							vValue[ip] = vValueOutside[ip]*vValueInside[ip]/(theta*vValueOutside[ip]+(1-theta)*vValueInside[ip]);
							//		                                UG_LOG("theta=" << theta << " vValue=" << vValue[ip] << "\n");
						}
					}
				}
			}
		}; // ip-cr-average
	}; // evaluate

	void update(){}

	  private:
	static const size_t max_number_of_ips = 20;
  };

template <typename TGridFunction>
class LevelSetUserVectorData
: public LevelSetUserDataBase<MathVector<TGridFunction::dim>, TGridFunction::dim,
  LevelSetUserVectorData<TGridFunction> >, virtual public INewtonUpdate
  {
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
	static const int dim = domain_type::dim;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	enum eval_type {sharp,cr_ip_average};

	  public:
	// set evaluation type, implemented so far:
	// 0 sharp (compute lsf in ip and give out inside value if inside or outside value if outside)
	// 1 cr_ip_average (compute averaged value in ip, if in CR-FV-Geometry from-value in ip and to-value in ip have different signs,
	//                  harmonic average is computed using intersection position from from and to node which is computed from level set function,
	//                  the ips given in evaluate must be the CR-FV ips)
	void set_eval_type(int type){
		if (type==0) m_eval_type=sharp;
		if (type==1) m_eval_type=cr_ip_average;
	}

	eval_type m_eval_type;

	  private:
	// level set grid function
	SmartPtr<TGridFunction> m_phi;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

	/*	//	component of function
			size_t m_fct;

		//	local finite element id
			LFEID m_lfeID;               */

	  public:
	void set_inside_data(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imInsideData = data;
	}

	void set_inside_data(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_inside_data(f);
	}

	void set_inside_data(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_inside_data(f);
		}
	}

	void set_inside_data(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_inside_data(f);
		}
	}

#ifdef UG_FOR_LUA
	void set_inside_data(const char* fctName)
	{
		set_inside_data(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
#endif

	void set_outside_data(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imOutsideData = data;
	}

	void set_outside_data(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_outside_data(f);
	}

	void set_outside_data(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_outside_data(f);
		}
	}

	void set_outside_data(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_outside_data(f);
		}
	}

#ifdef UG_FOR_LUA
	void set_outside_data(const char* fctName)
	{
		set_outside_data(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
#endif

	  private:
	///	Data import for inside and outside data
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imInsideData;
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imOutsideData;

	  public:
	/// constructor
	LevelSetUserVectorData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		m_phi = spGridFct;
		domain_type& domain = *m_phi->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		m_eval_type = cr_ip_average;
	}

	virtual ~LevelSetUserVectorData(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GridObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		const size_t numVertices = element->num_vertices();
		std::vector<Vertex*> childVertex(numVertices);
		std::vector<number> phi(numVertices);

		for (size_t i=0;i<nip;i++){
			//		            	UG_LOG("co(" << i << ",:)=" << vGlobIP[i] << "\n");
		}

		// find child vertices by injection
		for(size_t i = 0; i < numVertices; ++i){
			childVertex[i] = element->vertex(i);
		};

		// find out level, then get child vertices on finest level
		size_t numChildren = m_grid->template num_children<elem_type>(element);
		if (numChildren!=0){
			size_t lowerLevel=0;
			elem_type* childElem = m_grid->template get_child<elem_type>(element,0);
			numChildren = m_grid->template num_children<elem_type>(childElem);
			lowerLevel++;
			while (numChildren>0){
				childElem = m_grid->template get_child<elem_type>(childElem,0);
				numChildren = m_grid->template num_children<elem_type>(childElem);
				lowerLevel++;
			};
			for (size_t i=0;i<nip;i++){
				for (size_t j=0;j<lowerLevel;j++){
					childVertex[i]=m_grid->template get_child<Vertex>(childVertex[i],0);
				}
			}
		}
		//    create Multiindex
		std::vector<DoFIndex> ind;
		for (size_t i=0;i<numVertices;i++){
			m_phi->dof_indices(childVertex[i], 0, ind);
			phi[i]=DoFRef(*m_phi, ind[0]);
		};
		bool onls=false;
		bool inside=false;
		for (size_t i=0;i<numVertices;i++){
			if (phi[i]==0){
				continue;
			};
			if (phi[i]<0) inside=true;
			for (size_t j=i+1;j<numVertices;j++){
				if (phi[i]*phi[j]<0){
					onls = true;
					break;
				};
			}
		}
		if (onls==false){
			if (inside==false){
				(*m_imInsideData)(vValue,
						vGlobIP,
						time, si,
						elem,
						vCornerCoords,
						vLocIP,
						nip,
						u,
						vJT);
			} else {
				(*m_imOutsideData)(vValue,
						vGlobIP,
						time, si,
						elem,
						vCornerCoords,
						vLocIP,
						nip,
						u,
						vJT);
			};
			return;
		};
		MathVector<dim>  vValueInside[max_number_of_ips];
		MathVector<dim>  vValueOutside[max_number_of_ips];
		(*m_imInsideData)(vValueInside,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		(*m_imOutsideData)(vValueOutside,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		if (m_eval_type == sharp){
			for(size_t ip = 0; ip < nip; ++ip)
			{
				//	reference object id
				ReferenceObjectID roid = elem->reference_object_id();

				//	memory for shapes
				std::vector<number> vShape;

				// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
				const LocalShapeFunctionSet<refDim>& rTrialSpace =
						LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

				//  evaluate shapes at ip
				rTrialSpace.shapes(vShape, vLocIP[ip]);

				//    get multiindices of element
				std::vector<DoFIndex> ind;
				m_phi->dof_indices(elem, 0, ind);

				//     compute lsf at integration point
				number phiValue = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					phiValue += phi[sh] * vShape[sh];
				}
				if (phiValue<0) vValue[ip] = vValueInside[ip];
				else  vValue[ip] = vValueOutside[ip];
				//		                    UG_LOG("vValue(" << ip << ")=" << vValueOutside[ip] << "\n");
			};
		}

		if (m_eval_type == cr_ip_average){
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	memory for shapes
			std::vector<number> vShape;

			//    get domain of grid function
			const domain_type& domain = *m_phi->domain().get();

			//    get position accessor
			typedef typename domain_type::position_accessor_type position_accessor_type;
			const position_accessor_type& posAcc = domain.position_accessor();

			position_accessor_type aaPos = m_phi->domain()->position_accessor();

			//    coord and vertex array
			MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
			Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
			DimCRFVGeometry<dim> geo;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = element->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
			};
			//    evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			std::vector<number> phiSideValue(geo.num_sh());
			// compute interpolated lsf values in cr dofs
			for (size_t i=0;i<geo.num_scv();i++){
				// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
				const LocalShapeFunctionSet<dim>& rTrialSpace =
						LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

				//  evaluate shapes at ip
				rTrialSpace.shapes(vShape, geo.scv(i).local_ip());

				//    get multiindices of element
				std::vector<DoFIndex> ind;
				m_phi->dof_indices(elem, 0, ind);

				//     compute lsf at integration point
				phiSideValue[i] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					phiSideValue[i] += phi[sh] * vShape[sh];
				}
			}
			for (size_t ip=0;ip<nip;ip++){
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				number phiFrom = phiSideValue[scvf.from()];
				number phiTo = phiSideValue[scvf.to()];
				//		                    UG_LOG("ip=" << ip << " inside=" << vValueInside[ip] << " outside=" << vValueOutside[ip] << "\n");
				if (phiFrom==0){
					if (phiTo<=0){
						vValue[ip] = vValueInside[ip];
					} else {
						vValue[ip] = vValueOutside[ip];
					}
				} else {
					number theta;
					if (phiFrom<0){
						if (phiTo<=0){
							vValue[ip] = vValueInside[ip];
						} else {
							theta = phiFrom/(phiFrom-phiTo);
							for (int d=0;d<dim;d++)
								vValue[ip][d] = vValueInside[ip][d]*vValueOutside[ip][d]/(theta*vValueInside[ip][d]+(1-theta)*vValueOutside[ip][d]);
							//		                                UG_LOG("theta=" << theta << " vValue=" << vValue[ip] << "\n");
						}
					} else {
						if (phiTo>=0){
							vValue[ip] = vValueOutside[ip];
						} else {
							theta = phiFrom/(phiFrom-phiTo);
							for (int d=0;d<dim;d++)
								vValue[ip][d] = vValueOutside[ip][d]*vValueInside[ip][d]/(theta*vValueOutside[ip][d]+(1-theta)*vValueInside[ip][d]);
							//		                                UG_LOG("theta=" << theta << " vValue=" << vValue[ip] << "\n");
						}
					}
				}
			}
		}; // ip-cr-average
	}; // evaluate

	void update(){}

	  private:
	static const size_t max_number_of_ips = 20;
  };


template<typename TElem,typename TGrid>
void collect_finest_level_children(std::vector<TElem*>& childElemVector,TElem* elem,TGrid* grid){
	size_t numChildren = grid->template num_children<TElem>(elem);
	if (numChildren==0){
		childElemVector.push_back(elem);
	} else {
		for (size_t i=0;i<numChildren;i++){
			TElem* childElem = grid->template get_child<TElem>(elem,i);
			size_t nc=grid->template num_children<TElem>(childElem);
			if (nc==0){
				childElemVector.push_back(childElem);
			} else {
				collect_finest_level_children<TElem,TGrid>(childElemVector,childElem,grid);
			};
		};
	}
};

template<typename TElem,typename TPos,int dim>
void computeElemBarycenter(MathVector<dim>& bary,TElem* elem,TPos posA){
	bary = 0;
	size_t noc=elem->num_vertices();
	for (size_t i=0;i<noc;i++){
		bary+=posA[elem->vertex(i)];
	}
	bary/=noc;
}

template <typename TGridFunction>
class CRTwoPhaseSource
: public LevelSetUserDataBase<MathVector<TGridFunction::dim>, TGridFunction::dim,
  CRTwoPhaseSource<TGridFunction> >, virtual public INewtonUpdate
  {
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
	static const int dim = domain_type::dim;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	  private:
	// level set grid function
	SmartPtr<TGridFunction> m_phi;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

	  private:

	// gravitational constant
	number m_gravitational_constant;

	// surface tension factor
	number m_sigma;

	///	Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;

	  public:
	void set_gravitation(number gravityconst){
		m_gravitational_constant = gravityconst;
	}

	void set_sigma(number sigma){
		m_sigma = sigma;
	}

	/////////// Source

	void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imSource = data;
	}

	void set_source(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_source(f);
	}

	void set_source(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_source(f);
		}
	}

	void set_source(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_source(f);
		}
	}

#ifdef UG_FOR_LUA
	void set_source(const char* fctName)
	{
		set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
#endif

	  public:
	void set_density(SmartPtr<CplUserData<number, dim> > user){
		m_imDensity = user;
	}
	void set_density(number val){
		set_density(make_sp(new ConstUserNumber<dim>(val)));
	}
#ifdef UG_FOR_LUA
	void set_density(const char* fctName){
		set_density(LuaUserDataFactory<number, dim>::create(fctName));
	}
#endif

	  private:
	///	density import for inside and outside density
	SmartPtr<CplUserData<number,dim> > m_imDensity;

	  public:
	/// constructor
	CRTwoPhaseSource(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
		m_phi = spGridFct;
		for (int si=0;si<m_phi->num_subsets();++si){
			if (m_phi->num_fct(si)<2){
				UG_THROW("No curvature component in approximation space.");
			}
			if (m_phi->local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, dim,1)){
				UG_THROW("First component in approximation space must be of Lagrange 1 type.");
			}
			if (m_phi->local_finite_element_id(1) != LFEID(LFEID::PIECEWISE_CONSTANT,dim,0)){
				UG_THROW("Second component in approximation space must be of piecewise constant type.");
			}
		};
		domain_type& domain = *m_phi->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		set_gravitation(980);
		set_source(0.0);
		set_density(1);
	}

	virtual ~CRTwoPhaseSource(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GridObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		// evaluate source data
		//		   	   for (size_t i=0;i<nip;i++) vValue[i]=0;
		(*m_imSource)(vValue,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		// find corner values of level set function on finest level
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		const size_t numVertices = element->num_vertices();
		std::vector<Vertex*> childVertex(numVertices);
		std::vector<number> phi(numVertices);

		// find child vertices by injection
		for(size_t i = 0; i < numVertices; ++i){
			childVertex[i] = element->vertex(i);
		};

		bool onfinestlevel;
		// find out level, then get child vertices on finest level
		size_t numChildren = m_grid->template num_children<elem_type>(element);
		if (numChildren!=0){
			size_t lowerLevel=0;
			elem_type* childElem = m_grid->template get_child<elem_type>(element,0);
			numChildren = m_grid->template num_children<elem_type>(childElem);
			lowerLevel++;
			while (numChildren>0){
				childElem = m_grid->template get_child<elem_type>(childElem,0);
				numChildren = m_grid->template num_children<elem_type>(childElem);
				lowerLevel++;
			};
			for (size_t i=0;i<nip;i++){
				for (size_t j=0;j<lowerLevel;j++){
					childVertex[i]=m_grid->template get_child<Vertex>(childVertex[i],0);
				}
			}
		} else onfinestlevel=true;

		//    create Multiindex
		std::vector<DoFIndex> ind;
		for (size_t i=0;i<numVertices;i++){
			m_phi->dof_indices(childVertex[i], 0, ind);
			phi[i]=DoFRef(*m_phi, ind[0]);
		};
		bool onls=false;
		//bool inside=false; // inside seems to be unused (a.vogel)
		for (size_t i=0;i<numVertices;i++){
			//debug UG_LOG("phi(" << i+1 << ")=" << phi[i] << "\n");
			if (phi[i]==0){
				continue;
			};
			//if (phi[i]<0) inside=true;
			for (size_t j=i+1;j<numVertices;j++){
				//debug UG_LOG(" phi(" << j+1 << ")=" << phi[j] << "\n");
				if (phi[i]*phi[j]<0){
					onls = true;
					break;
				};
			}
			if (onls==true) break;
		}
		if (m_gravitational_constant!=0){
			number densityValue[max_number_of_ips];
			(*m_imDensity)(densityValue,
					vGlobIP,
					time, si,
					elem,
					vCornerCoords,
					vLocIP,
					nip,
					u,
					vJT);
			for (int d=0;d<dim;d++){
				for (size_t i=0;i<nip;i++){
					vValue[i][d] += m_gravitational_constant * densityValue[i];
				}
			}
		}
		if (onls==true){
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	memory for shapes
			std::vector<number> vShape;

			//    get domain of grid function
			const domain_type& domain = *m_phi->domain().get();

			//    get position accessor
			typedef typename domain_type::position_accessor_type position_accessor_type;
			const position_accessor_type& posAcc = domain.position_accessor();

			position_accessor_type aaPos = m_phi->domain()->position_accessor();

			//    coord and vertex array
			MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
			Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
			elem_type* evaluationElem;
			MathVector<dim> surftensSource[max_number_of_ips];
			for (size_t i=0;i<nip;i++){
				surftensSource[i]=0;
			}
			if (onfinestlevel==true){
				evaluationElem = element;
			} else {
				// compute barycenter
				MathVector<dim> coarseElementBary;
				computeElemBarycenter<elem_type,position_accessor_type,dim>(coarseElementBary,element,aaPos);
				// collect child elements
				std::vector<elem_type*> childElems;
				collect_finest_level_children<elem_type,grid_type>(childElems,element,m_grid);
				evaluationElem = childElems[0];
				number dist=1e+9;
				for (size_t ce=0;ce<childElems.size();ce++){
					MathVector<dim> bary;
					elem_type* currentElem = childElems[ce];
					size_t elnoc = currentElem->num_vertices();
					bool currentElemOnLS=false;
					std::vector<number> currentPhi(elnoc);
					for (size_t j=0;j<elnoc;j++){
						m_phi->dof_indices(currentElem->vertex(j),0,ind);
						currentPhi[j]=DoFRef(*m_phi,ind[0]);
					}
					for (size_t i=0;i<numVertices;i++){
						if (currentPhi[i]==0){
							continue;
						};
						for (size_t j=i+1;j<numVertices;j++){
							if (currentPhi[i]*currentPhi[j]<0){
								currentElemOnLS = true;
								break;
							}
						}
						if (currentElemOnLS==true) break;
					}
					if (currentElemOnLS==false) continue;
					computeElemBarycenter<elem_type,position_accessor_type,dim>(bary,currentElem,aaPos);
					number localdist  = VecDistance(bary,coarseElementBary);
					if (localdist<dist){
						dist = localdist;
						evaluationElem = childElems[ce];
					}
				}
			}
			m_phi->dof_indices(evaluationElem,1, ind);
			// evaluate curvature
			number element_curvature=DoFRef(*m_phi, ind[0]);
			// UG_LOG(element_curvature << "\n");
			// scale with surface tension factor
			element_curvature*=m_sigma;
			DimCRFVGeometry<dim> geo;
			const size_t eeNumVertices = evaluationElem->num_vertices();
			for(size_t i = 0; i < eeNumVertices; ++i){
				vVrt[i] = evaluationElem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				// UG_LOG("co(" << i+1 << ",:)=[" << coCoord[i][0] << "," << coCoord[i][1] << "];\n");
			};
			// evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			std::vector<number> phiSideValue(geo.num_sh());
			// compute interpolated lsf values in cr dofs
			for (size_t i=0;i<geo.num_scv();i++){
				// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
				const LocalShapeFunctionSet<dim>& rTrialSpace =
						LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim,1));

				//  evaluate shapes at ip
				rTrialSpace.shapes(vShape, geo.scv(i).local_ip());

				//    get multiindices of element
				std::vector<DoFIndex> ind;
				m_phi->dof_indices(elem, 0, ind);

				//     compute lsf at integration point
				phiSideValue[i] = 0.0;
				for(size_t sh = 0; sh < vShape.size(); ++sh)
				{
					phiSideValue[i] += phi[sh] * vShape[sh];
				}
			}
			//for (size_t i=0;i<geo.num_scv();i++){//debug
			//debug  UG_LOG("phi(" << i << ")=" << phiSideValue[i] << "\n");//debug
			//}//debug
			size_t nscvf = geo.num_scvf();
			for (size_t ip=0;ip<nscvf;ip++){
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				number phiFrom = phiSideValue[scvf.from()];
				number phiTo = phiSideValue[scvf.to()];
				if (phiFrom<0)   for (int d=0;d<dim;d++) surftensSource[scvf.from()][d] += scvf.normal()[d] * element_curvature;
				if (phiTo<0) for (int d=0;d<dim;d++) surftensSource[scvf.to()][d] -= scvf.normal()[d] * element_curvature;
			}
			for (size_t i=0;i<geo.num_scv();i++){
				surftensSource[i]/=geo.scv(i).volume();
				vValue[i] += surftensSource[i];
			}
			//for (size_t i=0;i<nip;i++){
				//	UG_LOG("rhs(" << i << ")=" << vValue[i] << "\n");//debug
			//};
			//UG_LOG("##############\n");//debug
		}
	}; // evaluate

	void update(){}

	  private:
	static const size_t max_number_of_ips = 20;
  };


}; // end namespace levelset
}; // end namespace ug

#endif /* __LEVEL_SET__LEVEL_SET_USER_DATA__ */
