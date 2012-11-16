/*
 * level_set_user_data.h
 *
 *  Created on: 15.11.2012
 *      Author: Christian Wehner
 */

#ifndef __LEVEL_SET__LEVEL_SET_USER_DATA__
#define __LEVEL_SET__LEVEL_SET_USER_DATA__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/lib_disc.h"


namespace ug{
namespace LevelSet{

template <typename TData, int dim, typename TImpl>
class LevelSetUserDataBase
	: 	public UserData<TData,dim>
{
	public:
		////////////////
		// one value
		////////////////
		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si) const
		{
			UG_THROW("LevelSetUserData: Need element.");
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<1>& locIP) const
		{
			getImpl().template evaluate<1>(&value,&globIP,time,si,u,elem,vCornerCoords,&locIP, 1, NULL);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<2>& locIP) const
		{
			getImpl().template evaluate<2>(&value,&globIP,time,si,u,elem,vCornerCoords,&locIP, 1, NULL);
		}

		virtual void operator() (TData& value,
		                         const MathVector<dim>& globIP,
		                         number time, int si,
		                         LocalVector& u,
		                         GeometricObject* elem,
		                         const MathVector<dim> vCornerCoords[],
		                         const MathVector<3>& locIP) const
		{
			getImpl().template evaluate<3>(&value,&globIP,time,si,u,elem,vCornerCoords,&locIP, 1, NULL);
		}

		////////////////
		// vector of values
		////////////////

		virtual void operator() (TData vValue[],
		                         const MathVector<dim> vGlobIP[],
		                         number time, int si, const size_t nip) const
		{
			UG_THROW("LevelSetUserData: Need element.");
		}


		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<1> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<1, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<1>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<2> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<2, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<2>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void operator()(TData vValue[],
		                        const MathVector<dim> vGlobIP[],
		                        number time, int si,
		                        LocalVector& u,
		                        GeometricObject* elem,
		                        const MathVector<dim> vCornerCoords[],
		                        const MathVector<3> vLocIP[],
		                        const size_t nip,
		                        const MathMatrix<3, dim>* vJT = NULL) const
		{
			getImpl().template evaluate<3>(vValue,vGlobIP,time,si,u,elem,
			                               vCornerCoords,vLocIP,nip, vJT);
		}

		virtual void compute(LocalVector* u, GeometricObject* elem, bool bDeriv = false)
		{
			const number t = this->time();
			const int si = this->subset();
			for(size_t s = 0; s < this->num_series(); ++s)
				getImpl().template evaluate<dim>(this->values(s), this->ips(s), t, si,
												*u, elem, NULL, this->template local_ips<dim>(s),
												this->num_ip(s));
		}

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
		typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;

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
	/// constructor
		LevelSetUserData(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct){
			m_phi = spGridFct;
			domain_type& domain = *m_phi->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
			m_spApproxSpace = approxSpace;
		}

		virtual ~LevelSetUserData(){};

			template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
			elem_type* element = static_cast<elem_type*>(elem);

			const size_t numVertices = element->num_vertices();
			std::vector<VertexBase*> childVertex(numVertices);
			std::vector<number> phi(numVertices);

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
						childVertex[i]=m_grid->template get_child<VertexBase>(childVertex[i],0);
					}
				}
			}
			//	create Multiindex
			std::vector<MultiIndex<2> > ind;
			for (size_t i=0;i<numVertices;i++){
				m_phi->multi_indices(childVertex[i], 0, ind);
				phi[i]=DoFRef(*m_phi, ind[0]);
			};
				/*	bool onls=false;
			bool inside=false;
			for (size_t i=0;i<numVertices;i++){
				if (phi[i]==0){
			        continue;
				};
				if (phi[i]<0) inside=true;
			    for (size_t j=i+1;j<elementnoc;j++){
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
			                    u,
			                    elem,
			                    vCornerCoords,
			                    vLocIP,
			                    nip,
			                    vJT);
				} else {
					(*m_imOutsideData)(vValue,
			                    vGlobIP,
			                    time, si,
			                    u,
			                    elem,
			                    vCornerCoords,
			                    vLocIP,
			                    nip,
			                    vJT);
				};
				return;
			};
			(*m_imInsideData)(vValueInside,
						                    vGlobIP,
						                    time, si,
						                    u,
						                    elem,
						                    vCornerCoords,
						                    vLocIP,
						                    nip,
						                    vJT);
			(*m_imOutsideData)(vValueOutside,
						                    vGlobIP,
						                    time, si,
						                    u,
						                    elem,
						                    vCornerCoords,
						                    vLocIP,
						                    nip,
						                    vJT);
							if (m_eval_type == sharp){
				for(size_t ip = 0; ip < nip; ++ip)
				{
					// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
					const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalShapeFunctionSetProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, 1));

					//	evaluate at shapes at ip
					rTrialSpace.shapes(vShape, vLocIP[ip]);

					//	get multiindices of element
					std::vector<MultiIndex<2> > ind;
					m_spGridFct->multi_indices(elem, 0, ind);

					// 	compute lsf at integration point
					phiValue = 0.0;
					for(size_t sh = 0; sh < vShape.size(); ++sh)
					{
						const number valSH = DoFRef( *m_spGridFct, ind[sh]);
						phiValue += phi[sh] * vShape[sh];
					}
					if (phiValue<0) vValue[ip] = vInside[ip];
						else  vValue[ip] = vOutside[ip];
				};
			}
			if (m_eval_type == cr_ip_average){
				//	get domain of grid function
				domain_type& domain = *m_u->domain().get();

				//	get position accessor
				typedef typename domain_type::position_accessor_type position_accessor_type;
				const position_accessor_type& posAcc = domain.position_accessor();

				position_accessor_type aaPos = m_phi->domain()->position_accessor();

				//	coord and vertex array
				MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
				VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
				DimCRFVGeometry<dim> geo;

				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
					// UG_LOG("co_coord(" << i<< "+1,:)=" << coCoord[i] << "\n");
				};
				//	evaluate finite volume geometry
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				// compute lsf values in cr dofs
				for (size_t i=0;i<geo.num_scv();i++){
					// compute lsf value in ip (from Lagrange-1 shape function) and set value according to position (inside or outside)
					const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalShapeFunctionSetProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, 1));

					//	evaluate at shapes at ip
					rTrialSpace.shapes(vShape, geo->scv(i).local_ip());

					//	get multiindices of element
					std::vector<MultiIndex<2> > ind;
					m_spGridFct->multi_indices(elem, 0, ind);

					// 	compute lsf at integration point
					phiSideValue[i] = 0.0;
					for(size_t sh = 0; sh < vShape.size(); ++sh)
					{
						const number valSH = DoFRef( *m_spGridFct, ind[sh]);
						phiSideValue[i] += phi[sh] * vShape[sh];
					}
				}
				for (size_t ip=0;ip<nip;ip++){
					phiFrom = phiSideValue[scvf.from()];
					phiTo = phiSideValue[scvf.to()];
					if (phiFrom==0){
						if (phiTo<=0){
							vValue[ip] = vValueInside[ip];
						} else {
							vValue[ip] = vValueOutside;
						}
					} else {
						if (phiFrom<0){
							if (phiTo<=0){
								vValue[ip] = vValueInside[ip];
							} else {
								theta = phiFrom/(phiFrom-phiTo);
								vValue[ip] = vValueInside[ip]*vValueOutside[ip]/(theta*vValueInside[ip]+(1-theta)*vValueOutside[ip]);
							}
						} else {
							if (phiTo>=0){
								vValue[ip] = vValueOutside[ip];
							} else {
								theta = phiFrom/(phiFrom-phiTo);
								vValue[ip] = vValueOutside[ip]*vValueInside[ip]/(theta*vValueOutside[ip]+(1-theta)*vValueInside[ip]);
							}
						}
					}
				}
		}; // ip-cr-average  */
	}; // evaluate

	void update(){}
};


}; // end namespace levelset
}; // end namespace ug

#endif /* __LEVEL_SET__LEVEL_SET_USER_DATA__ */
