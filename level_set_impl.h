/*
 * level_set_impl.h
 *
 *  Created on: 01.07.2011
 *      Author: Christian Wehner  christian.wehner@gcsc.uni-frankfurt.de
 */

#ifndef LEVEL_SET_UTIL_IMPL_H_
#define LEVEL_SET_UTIL_IMPL_H_

#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_grid/algorithms/attachment_util.h"

#include <algorithm>

namespace ug{
namespace LevelSet{

/**
 * limit previously computed gradient so that the control-volume-wise linear
 * interpolation function does not introduce new maxima or minima into the data 
 */
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::limit_grad
(
	TGridFunction& uOld,
	aaGrad& aaGrad
)
{
	// get grid
	grid_type& grid = *uOld.domain()->grid();

	position_accessor_type& aaPos = uOld.domain()->position_accessor();

	std::vector<DoFIndex> ind;

	//	create Attachment for scv-volume size
	ANumber aMax;
	ANumber aMin;

	//	attach to grid
	grid.attach_to_vertices(aMin);
	grid.attach_to_vertices(aMax);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaMin(grid, aMin);
	Grid::VertexAttachmentAccessor<ANumber> aaMax(grid, aMax);
	for (int si=0;si<uOld.num_subsets();++si)
	{
		for(VertexConstIterator iter = uOld.template begin<Vertex>(si);
				iter !=uOld.template end<Vertex>(si); ++iter)
		{
			Vertex* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
			//	read indices on vertex
			//	get vector holding all indices on the vertex
			uOld.inner_dof_indices(vrt, 0, ind);
			aaMax[vrt] = DoFRef(uOld, ind[0]);
			aaMin[vrt] = DoFRef(uOld, ind[0]);
		}
	}
	for (int si=0;si<uOld.num_subsets();++si)
	{
		//UG_LOG("si " << si << "\n");
		for(EdgeConstIterator iter = uOld.template begin<Edge>(si) ;
				iter !=uOld.template end<Edge>(si); ++iter)
		{
			Edge* edge = *iter;
			Vertex* vi=edge->vertex(0);
			Vertex* vj=edge->vertex(1);
			uOld.inner_dof_indices(vi, 0, ind);
			number ui = DoFRef(uOld, ind[0]);
			uOld.inner_dof_indices(vj, 0, ind);
			number uj = DoFRef(uOld, ind[0]);
			//UG_LOG("edge " << aaPos[vi] << "-" << aaPos[vj] << " [" << ui << " " << uj << "]\n");
			if (uj<aaMin[vi])
				aaMin[vi]=uj;
			if (uj>aaMax[vi])
				aaMax[vi]=uj;
			if (ui<aaMin[vj])
				aaMin[vj]=ui;
			if (ui>aaMax[vj])
				aaMax[vj]=ui;
		}

	};
/*
	for (int si=0;si<2;++si)
	{
		for(VertexConstIterator iter = uOld.template begin<Vertex>(si) ;
				iter !=uOld.template end<Vertex>(si); ++iter)
		{
		    Vertex* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
			// read indices on vertex
			//	get vector holding all indices on the vertex
			uOld.inner_dof_indices(vrt, 0, ind);
			//UG_LOG(coord << " min=" << aaMin[vrt] << " max=" << aaMax[vrt] << "\n");
		}
	}
*/
	for (int si=0;si<1;++si)
	{
		for(EdgeConstIterator iter = uOld.template begin<Edge>(si) ;
				iter !=uOld.template end<Edge>(si); ++iter)
		{
			Edge* edge = *iter;
			Vertex* vi=edge->vertex(0);
			Vertex* vj=edge->vertex(1);
			MathVector<dim> coordi,coordj,coordij,distVec,gradi,gradj;
			gradi = aaGrad[vi];
			gradj = aaGrad[vj];
			coordi = aaPos[vi];
			coordj = aaPos[vj];
			uOld.inner_dof_indices(vi, 0, ind);
			number ui = DoFRef(uOld, ind[0]);
			uOld.inner_dof_indices(vj, 0, ind);
			number uj = DoFRef(uOld, ind[0]);
			VecScaleAdd(coordij,0.5,coordi,0.5,coordj);
			VecSubtract(distVec, coordij,coordi);
			number uij = ui + distVec*gradi;
			number alpha = 1;
			if (uij>ui)
			{
				if (uij>aaMax[vi]) alpha=(aaMax[vi]-ui)/(distVec*gradi);
				if (alpha<1)
				{
					//UG_LOG("edge " << coordi << " " << coordj << "\n");
					//UG_LOG(coordi << " u " << ui << " uij "  << uij << " max " << aaMax[vi] << " alpha " << alpha << "\n");
					aaGrad[vi]*=alpha;
				};
			}
			else
			{
				if (uij<aaMin[vi]) alpha=(aaMin[vi]-ui)/(distVec*gradi);
				if (alpha<1)
				{
					//UG_LOG("edge " << coordi << " " << coordj << "\n");
					//UG_LOG(coordi << " u " << ui << " uij "  << uij << " min " << aaMax[vi] << " alpha " << alpha << "\n");
					aaGrad[vi]*=alpha;
				};
			};
			VecSubtract(distVec, coordij,coordj);
			uij = uj + distVec*gradj;
			alpha = 1;
			if (uij>uj)
			{
				if (uij>aaMax[vj]) alpha=(aaMax[vj]-uj)/(distVec*gradj);
				if (alpha<1)
				{
					//UG_LOG("-- edge " << coordi << " " << coordj << "\n");
					//UG_LOG(coordj << " u " << uj << " uij "  << uij << " max " << aaMax[vj] << " alpha " << alpha << "\n");
					aaGrad[vj]*=alpha;
				};
			}
			else
			{
				if (uij<aaMin[vj]) alpha=(aaMin[vj]-uj)/(distVec*gradj);
				if (alpha<1)
				{
					//UG_LOG("-- edge " << coordi << " " << coordj << "\n");
					//UG_LOG(coordj << " u " << uj << " uij "  << uij << " min " << aaMax[vj] << " alpha " << alpha << "\n");
					aaGrad[vj]*=alpha;
				};
			};
         // UG_LOG(" coord vertex 0 " << aaPos[v0] << " coord vertex 1 " << aaPos[v1] << "\n");
		}

		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
		//	coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		//	coord and vertex array
		MathVector<dim> grad[domain_traits<dim>::MaxNumVerticesOfElem];
		//	values of the grid function
		number u[domain_traits<dim>::MaxNumVerticesOfElem];
		//	get iterators
		ElemIterator iter = uOld.template begin<ElemType>(si);
		ElemIterator iterEnd = uOld.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;

			//	get position accessor
			const position_accessor_type& aaPos = uOld.domain()->position_accessor();

			//	compute center of mass
			MathVector<dim> center;
			center=0;
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
			{
				vVrt[i] = elem->vertex(i);
				coCoord[i] = aaPos[vVrt[i]];
				grad[i] = aaGrad[vVrt[i]] ;
				VecAppend(center,coCoord[i]);
				uOld.inner_dof_indices(vVrt[i], 0, ind);
				u[i]=DoFRef(uOld, ind[0]);
			};
			center/=noc;
			for (size_t i=0;i<noc;++i)
			{
				number alpha=1;
				MathVector<dim> distVec;
				number uCenter;
				VecSubtract(distVec,center,coCoord[i]);
				uCenter = u[i] + distVec*grad[i];
				if (uCenter>u[i])
				{
					if (uCenter>aaMax[vVrt[i]]) alpha=(aaMax[vVrt[i]]-u[i])/(distVec*grad[i]);
					if (alpha<1)
					{
						// UG_LOG("* " << coCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " max " << aaMax[vVrt[i]] << " alpha " << alpha << "\n");
						aaGrad[vVrt[i]]*=alpha;
					};
				}
				else
				{
					if (uCenter<aaMin[vVrt[i]]) alpha=(aaMin[vVrt[i]]-u[i])/(distVec*grad[i]);
					if (alpha<1)
					{
						// UG_LOG("*#* " << coCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " min " << aaMin[vVrt[i]] << " alpha " << alpha << "\n");
						aaGrad[vVrt[i]]*=alpha;
					};
				};
			};
		};
	};
	//	detach from grid
	grid.detach_from_vertices(aMin);
	grid.detach_from_vertices(aMax);
	return true;
};


/**
 * assemble element using upwind
 *
 * This function computes the contribution of the explicit local discretization
 * in one grid element.
 */
template<typename TGridFunction>
template <typename TElem>
bool FV1LevelSetDisc<TGridFunction>::assemble_element
(
	TElem& elem, ///< the element to compute the contribution for
	DimFV1Geometry<dim>& geo,
	grid_type& grid,
	TGridFunction& uNew,
	const TGridFunction& uOld,
	aaGrad& aaGradient,
	aaVol& aaVolume
)
{
//	a large enough number for max. number of corners and subcontrol volume faces
	static const size_t maxNumCo = 20;
	
//	get domain
	domain_type& domain = *uNew.domain().get();

//	create Multiindex
	std::vector<DoFIndex> multInd;

	//	hard code function (fct=0)
	//\todo: generalize
//	size_t fct=0;

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor();

//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

//	get vertices and extract corner coordinates
	const size_t numVertices = elem->num_vertices();
	for(size_t i = 0; i < numVertices; ++i)
	{
		vVrt[i] = elem->vertex(i);
		coCoord[i] = aaPos[vVrt[i]];
	};
	
// update fv geometry
	geo.update(elem, &(coCoord[0]), uOld.domain()->subset_handler().get());

//  fill node value vector
	std::vector<number> uValue(geo.num_scv());
	size_t noc = geo.num_scv();
	for (size_t i=0; i < noc; i++)
	{
		// if (dd.template inner_dof_indices<Vertex>(vVrt[i], 0, multInd) != 1) return false;
		uOld.inner_dof_indices(vVrt[i], 0, multInd);
		uValue[i]=DoFRef(uOld, multInd[0]);
	}
//  fill grad vector
    MathVector<dim> grad[maxNumCo];
	for (size_t i=0; i < noc; i++)
	{
		// if (dd.template inner_dof_indices<Vertex>(vVrt[i], 0, multInd) != 1) return false;
		uOld.inner_dof_indices(vVrt[i], 0, multInd);
		grad[i]= aaGradient[vVrt[i]];
	}
	
//  fill corner velocity vector
	MathVector<dim> coVelocity[maxNumCo];
	for (size_t i=0; i<noc; ++i)
		coVelocity[i]=0;
	if (m_gamma!=0)
	{
		const int si = 0; //TODO this should be corrected
		if (m_imVelocity->requires_grid_fct())
		{
			//	create storage
			LocalIndices localind;
			LocalVector localu;

			// 	get global indices
			uNew.indices(elem, localind);

			// 	adapt local algebra
			localu.resize(localind);

			// 	read local values of u
			GetLocalVector(localu, uNew);

			(*m_imVelocity)(coVelocity, geo.scv_global_ips(), m_time, si,
				elem,coCoord, geo.scv_local_ips(),geo.num_sh(), &localu);
		}
		else
		{
			// see user_data.h : 410
			(*m_imVelocity)(coVelocity, geo.scv_global_ips(), m_time, si, geo.num_sh());
		}
        if (m_gamma!=1) for (size_t i=0; i < noc; i++) coVelocity[i]*=m_gamma;
	}
    if (m_delta!=0)
    {
    	for (size_t i=0; i < noc; i++)
    	{
    	    number vnorm = VecLength(grad[i]);
    	    if (vnorm>1e-15) for (int j=0; j < dim; j++) coVelocity[i][j] += m_delta/vnorm*grad[i][j];
    	};
    };
    
//  compute the ip velocity from the corner velocity
	MathVector<dim> ipVelocity[maxNumCo];
	for (size_t ip=0;ip < geo.num_scvf();ip++)
	{
		ipVelocity[ip] = 0;
		const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
		for (size_t co=0; co < noc; co++)
		{
			for (int j=0; j<dim; j++)
			    ipVelocity[ip][j] += scvf.shape(co)*coVelocity[co][j];
		};
	}
	
//  fill source vector
	std::vector<number> coSource(noc);
	const int si = 0; //TODO this should be corrected
	(*m_imSource)(&coSource[0], geo.scv_global_ips(), m_time, si, geo.num_sh());

// compute fluxes
	size_t base;
	number flux;
	MathVector<dim> distVec;
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	    MathVector<dim> bNode;
	    // 	get current SCVF
	    const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
	    MathVector<dim>	ipCoord = scvf.global_ip();
	    size_t from = scvf.from();
	    size_t to   = scvf.to();
	    if (scvf.normal()*ipVelocity[ip]>0)
		    base = from;
		else
		    base = to;
		VecSubtract(distVec, ipCoord,coCoord[base]);
		// flux = v * n * u_{ip(i)}^{n+0.5}
		flux = m_dt*(ipVelocity[ip]*scvf.normal())*( uValue[base] + (distVec*grad[base]) + 0.5*m_dt*(coSource[base] - (grad[base]*coVelocity[base])) );
		//UG_LOG(ip << " flux=" << flux << "\n");
		uOld.inner_dof_indices(vVrt[from], 0, multInd);
        DoFRef(uNew, multInd[0])-=flux/aaVolume[ vVrt[from] ];
		if (! m_divFree)
		{
		    DoFRef(uNew, multInd[0])+=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ];
		    //DofRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.0*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ];
		};
		uOld.inner_dof_indices(vVrt[to], 0, multInd);
        DoFRef(uNew, multInd[0])+=flux/aaVolume[ vVrt[to] ];
		if (! m_divFree)
		{
		    DoFRef(uNew, multInd[0])-=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to])))/aaVolume[ vVrt[to] ];
		};
        number localCFL = std::max(m_dt*std::abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[from] ],m_dt*std::abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[to] ] );
        if (localCFL>m_maxCFL)
        {
            m_maxCFL = localCFL;
        };
	};
	
// boundary
	if (geo.num_bf()>0)
	{
		for (int si=0;si<uNew.num_subsets();++si)
		{
			//UG_LOG("si=" << si << "  num_bf(si)=" << geo.num_bf(si) << "\n");
			for(size_t i = 0; i < geo.num_bf(si); ++i)
			{
			// 	get current BF
				const typename DimFV1Geometry<dim>::BF& bf = geo.bf(si, i);
				const size_t nodeID = bf.node_id();
				MathVector<dim> bipVelocity;
				number bipU=0;
				number bipSource=0;
				MathVector<dim> bipGrad;
				bipVelocity=0;
				bipGrad=0;
				const MathVector<dim>* globalGradVec = bf.global_grad_vector();
				bipVelocity=0;
				for (size_t co=0;co<noc;co++)
				{
					//UG_LOG("corner " << co << " num_sh " << bf.num_sh() << "\n");
					for (int j=0;j<dim;j++)
					{
						bipVelocity[j] += bf.shape(co) * coVelocity[co][j];
						bipGrad[j] += bf.shape(co) * globalGradVec[co][j];
					};
					bipU += bf.shape(co) * uValue[co];
					bipSource += bf.shape(co) * coSource[co];
				}
				//flux = m_dt*(bipVelocity*bf.normal())*( bipU + 0.5*m_dt*(bipSource - (bipGrad*bipVelocity)) );
				flux = m_dt*(bipVelocity*bf.normal())*uValue[nodeID];// first order approximation
				uOld.inner_dof_indices(vVrt[nodeID], 0, multInd);
				DoFRef(uNew, multInd[0])-=flux/aaVolume[ vVrt[nodeID] ];
				if (!m_divFree)
				{
					DoFRef(uNew, multInd[0])+=m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] )/aaVolume[ vVrt[nodeID] ];// first order approximation
				};
			};
		};
	};
	
	return true;
}

/**
 * computes gradient in grid vertices as well as volume of control volumes
 * and saves the result in the grid attachments
 */
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::calculate_vertex_grad_vol
(
	TGridFunction& u,
	aaGrad& aaGradient,
	aaVol& aaVolume
)
{
//	get domain
	domain_type& domain = *u.domain().get();

//	get grid of domain
	typename domain_type::grid_type& grid = *domain.grid();

//	create Multiindex
	std::vector<DoFIndex> multInd;

//	create a FV Geometry for the dimension
	DimFV1Geometry<dim> geo;

//	hard code function (fct=0)
//\todo: generalize
	size_t fct=0;

	// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

//	sum up all contributions of the sub control volumes to one vertex in an attachment
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
	//	skip boundary
		if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)) continue;
		if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)) continue;
		if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)) continue;
		
	//	get iterators
		ElemIterator iter = u.template begin<ElemType>(si);
		ElemIterator iterEnd = u.template end<ElemType>(si);

	//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
		//	get Elem
			ElemType* elem = *iter;

		//	get position accessor
			const position_accessor_type& aaPos = domain.position_accessor();

		//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i)
			{
				vVrt[i] = elem->vertex(i);
				coCoord[i] = aaPos[vVrt[i]];
			};

		//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			//UG_LOG("Num Verts loaded: "<<vVrt.size()<<"\n");
			//UG_LOG("Num SCV computed: "<<geo.num_scv()<<"\n");

			number uValue[domain_traits<dim>::MaxNumVerticesOfElem];

		//	read indices on vertex
			size_t noc = geo.num_scv();
			for (size_t i=0;i < noc;i++)
			{
			//	get indices of function fct on vertex
				u.inner_dof_indices(vVrt[i], fct, multInd);

			//	read value of index from vector
				uValue[i]=DoFRef(u, multInd[0]);

			//	debug log
				//UG_LOG("corner " << i << " " << uValue[i] << "\n");
			}

		//	storage for global gradient
			MathVector<dim> globalGrad;

		//	loop corners
			for (size_t i=0;i < noc;i++)
			{
			//	get scv for sh
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(i);

			//	debug log
				//UG_LOG("gradient for corner " << i << "\n");

			//	reset global gradient
				globalGrad = 0.0;

			//	sum up gradients of shape functions in corner
				for(size_t sh = 0 ; sh < noc; ++sh)
				{
					//UG_LOG("local grad " << sh << " : " << scv.local_grad(sh) << "\n");
					//UG_LOG("unscaled global grad " << sh << " = " << scv.global_grad(sh) << "\n");
					//UG_LOG("uvalue(" << sh << ") =" << uValue[sh] << "\n");

					VecScaleAppend(globalGrad, uValue[sh], scv.global_grad(sh));
				}

			//	volume of scv
				number vol = scv.volume();

				//UG_LOG("*** global grad " << i << ": " << globalGrad << "\n");

			//	scale gradient by volume
				globalGrad *= vol;

			//	add both values to attachements
				aaGradient[vVrt[i]] += globalGrad;
				aaVolume[vVrt[i]] += vol;
			};
		}
	}

//	divide the gradients by the volumes
	// int count=0;
	position_accessor_type aaPos = u.domain()->position_accessor();
	for (int si=0;si < u.num_subsets();++si)
	{
		//UG_LOG("si " << si << "\n");
	    for(VertexConstIterator iter = u.template begin<Vertex>(si);
						   iter != u.template end<Vertex>(si); ++iter)
	    {
	    //	get vertex
		    Vertex* vrt = *iter;
		    if (aaVolume[vrt]!=0)
		    {
		        (aaGradient[vrt]) /= aaVolume[vrt];
		    };
		    //exact[0] = cos(coord[0]);// 6*coord[0];
		    //exact[1] = -4*sin(coord[1]);//-4*coord[1];
		    //number gError = sqrt( (exact[0]-aaGradient[vrt][0])*(exact[0]-aaGradient[vrt][0]) + (exact[1]-aaGradient[vrt][1])*(exact[1]-aaGradient[vrt][1]) );
	        //UG_LOG(count << "[ " << coord[0] << "," << coord[1] << " ] vol= " << aaVolume[vrt] << " " << "grad= ["
	        //		<< aaGradient[vrt][0] << "," << aaGradient[vrt][1] << "] exact grad = ["
	        //		<< exact[0] << "," << exact[1] << "] error: " << gError <<  "\n");
	        // count++;
	    }
	}
	
	//UG_LOG("#*#*#*#\n");
	return true ;
}

/**
 * sets Dirichlet values in solution vector for vertices in a given subset
 */
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::assign_dirichlet
(
	TGridFunction& numsol
)
{
//	get domain of grid function
	domain_type& domain = *numsol.domain().get();

//	UG_LOG("dirichlet\n");

	// UG_LOG("nr dir ss " << m_dirichlet_sg.size() << "\n");

	for(size_t i = 0; i < m_dirichlet_sg.size(); ++i)
	{
		const int si = m_dirichlet_sg[i];
		// UG_LOG("Dirichlet boundary is: "<<si<< "\n");
		for(VertexConstIterator iter = numsol.template begin<Vertex>(si);
									   iter != numsol.template end<Vertex>(si); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
			number exactVal;
			position_accessor_type aaPos = domain.position_accessor();

		//	get vector holding all indices on the vertex
			std::vector<DoFIndex> ind;

			const size_t numInd = numsol.inner_dof_indices(vrt, 0, ind);

		//	check indices
			if(numInd != 1) {UG_LOG("ERROR: Wrong number of indices!"); return false;}

			(*m_imDirichlet)(&exactVal,&aaPos[vrt],m_time,si,1);
			DoFRef(numsol, ind[0]) = exactVal;

			// MathVector<dim> coord = aaPos[vrt];
			// UG_LOG("coord " << coord[0] << "," << coord[1] << " <> " << BlockRef(numsol[ind[0][0]],ind[0][1]) << "\n");

			//if ((coord[0]==-1)||(coord[0]==1)||(coord[1]==-1)||(coord[1]==1)){
				//BlockRef(numsol[ind[0][0]],ind[0][1]) = exactVal;
			//};
			//if ((coord[0]==0)||(coord[0]==1)||(coord[1]==0)||(coord[1]==1)){
			//};
		 }
	}
	
	return true;
}

/**
 * Computes the time steps of the discretization of the level-set equation.
 *
 * Main function for assembling and solving ls equation as described in 
 * Frolkovic/Mikula, HIGH-RESOLUTION FLUX-BASED LEVEL SET METHOD, SIAM
 */
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::advect_lsf
(
	TGridFunction& uNew,
	TGridFunction& uOld
)
{
//	get domain of grid function
	domain_type& domain = *uNew.domain().get();

//	get grid of domain
	grid_type& grid = *domain.grid();

//	create Attachment for scv-volume size
	ANumber aScvVolume;

//	typedef of gradient attachment
	typedef Attachment<MathVector<dim> > AGradient;

//	create Attachment for gradient
	AGradient aGradient;

//	attach to grid
	grid.attach_to_vertices(aScvVolume);
	grid.attach_to_vertices(aGradient);

//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaVolume(grid, aScvVolume);
	Grid::VertexAttachmentAccessor<AGradient> aaGradient(grid, aGradient);

// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);


//	the CFL constant to compute
	m_maxCFL=0;
	
	//UG_LOG("***************************************************\n");
    //UG_LOG("***************************************************\n");
	//UG_LOG("***************************************************\n");

	VecAssign(uNew,uOld);
    MathVector<dim> coord;
	std::vector<DoFIndex> ind;
    position_accessor_type aaPos = domain.position_accessor();

//	compute time steps
	for (size_t step=0; step<m_nrOfSteps; step++)
	{
	// calculate scv volume and the gradient
	    if (! calculate_vertex_grad_vol(uNew,aaGradient, aaVolume)) {UG_LOG("ERROR: gradient computation failed!"); };
	    // SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0); // for debug set gradient to 0
	    if (m_limiter)
	    	limit_grad(uNew,aaGradient);
	    
	    // UG_LOG("num_subsets: " << uOld.num_subsets() << "\n");
	    
	//	loop over subsets to compute the new solution
	    for (int si=0;si<uOld.num_subsets();++si)
	    {
	    //	skip boundaries
	        //UG_LOG("si " << si << "\n");
	        if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)) continue;
	        if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)) continue;
	        if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)) continue;
	        
	        //UG_LOG("... \n");
		    //UG_LOG("***************************************************\n");
	        //UG_LOG("***********************" << si << "**************************\n");
		    //UG_LOG("***************************************************\n");
		    
		//	get iterators
		    ElemIterator iter = uNew.template begin<ElemType>(si);
		    ElemIterator iterEnd = uNew.template end<ElemType>(si);

		    DimFV1Geometry<dim> geo;

		// flag given neumann bnd subsets at the geometry, such that the
		//	geometry produces boundaryfaces (BF) for all sides of the
		//	element, that is in one of the subsets
		    for(size_t i = 0; i < m_neumann_sg.size(); ++i)
		    {
		        const int bndSi = m_neumann_sg[i];
//		        UG_LOG("Neumann boundary is: "<<bndSi<< "\n");
		        geo.add_boundary_subset(bndSi);
		    }

		//	loop elements compute the new solution
		    for(  ;iter !=iterEnd; ++iter)
		    {
		        //	get Elem
			    //UG_LOG("*** ELEM ***\n");
			    ElemType* elem = *iter;
			    //UG_LOG("element \n");
			    // uNew = uOld
			    assemble_element(elem, geo, grid, uNew, uOld, aaGradient, aaVolume);
		    };
	    };
	    
	//	take into account the source at the vertices
		for (int si=0;si<uOld.num_subsets();++si)
		{
			if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)) continue;
			if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)) continue;
			std::vector<number> sourceValue(1);
			MathVector<dim> sourceCo[1];
			for(VertexConstIterator iter = uNew.template begin<Vertex>(si);
					    			                      iter != uNew.template end<Vertex>(si); ++iter)
			{
				Vertex* vrt = *iter;
				sourceCo[0]= aaPos[vrt];
				(*m_imSource)(&sourceValue[0],sourceCo,m_time,si,1);
				uNew.inner_dof_indices(vrt, 0, ind);
				DoFRef(uNew, ind[0]) += m_dt*sourceValue[0];
			};
	    }
	    
	//	the new solution computed
	    m_time += m_dt;
	    m_timestep_nr++;
	    
	//	set the Dirichlet values
	    assign_dirichlet(uNew);
	    
	//	overwrite inactive nodes with old solution
	    for(size_t i = 0; i < m_inactive_sg.size(); ++i)
	    {
	    	 const int si = m_inactive_sg[i];
	    	 UG_LOG("inactive si: " << si << "\n");
	    	 for(VertexConstIterator iter = uNew.template begin<Vertex>(si);
	    	 									   iter != uNew.template end<Vertex>(si); ++iter)
	    	 {
	    	     Vertex* vrt = *iter;
	    	     UG_LOG("*\n");
	    	     uNew.inner_dof_indices(vrt, 0, ind);
	    		 DoFRef(uNew, ind[0]) = DoFRef(uOld, ind[0]);
	    	 }
	    };
	
	//	time step done:
	    UG_LOG("time step length: " << m_dt << "\n");
	    UG_LOG("time step no.: " << m_timestep_nr << "\n");
	    UG_LOG("time: " << m_time << "\n");
        UG_LOG("max CFL: " << m_maxCFL << "\n");
        
	    if (m_nrOfSteps>1)
	    {
	    	// Attention, uOld is overwritten
	    	VecAssign(uOld,uNew);
	    };
	};
	
    //	detach from grid
	grid.detach_from_vertices(aScvVolume);
	grid.detach_from_vertices(aGradient);

//	done
	return true;
}

/**
 * compute volume of control volumes and saves it in the attachment
 */
template <typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::calculate_vertex_vol
(
	TGridFunction& u,
	aaVol& aaVolume
)
{
//	get domain
	domain_type& domain = *u.domain().get();

//	create a FV Geometry for the dimension
	DimFV1Geometry<dim> geo;

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor();

//	sum up all contributions of the sub control volumes to one vertex in an attachment
	for(int si = 0; si < u.num_subsets(); ++si)
	{
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	//	get iterators
		ElemIterator iter = u.template begin<ElemType>(si);
		ElemIterator iterEnd = u.template end<ElemType>(si);

	//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
		//	get Elem
			ElemType* elem = *iter;

			size_t numVertices=elem->num_vertices();
		//	extract corner coordinates
			for(size_t i = 0; i < numVertices; ++i)
				coCoord[i] = aaPos[elem->vertex(i)];

		//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());


		//	loop corners
			for (size_t i=0;i < geo.num_scv();i++)
			{
				//	get scv for sh
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(i);

				aaVolume[elem->vertex(i)] += scv.volume();
			}
		}
	}

	return true;
}

/**
 * computes error w.r.t. the analytical solution (taken from the dirichlet bc)
 *
 * This is used in test cases.
 */
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_error(TGridFunction& numsol)
{
//	get domain of grid function
	domain_type& domain = *numsol.domain().get();

//	get grid of domain
	grid_type& grid = *domain.grid();

	//	create Attachment for scv-volume size
	ANumber aScvVolume;

	//	attach to grid
	grid.attach_to_vertices(aScvVolume);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaVolume(grid, aScvVolume);

	// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);

    number l1Error=0;
	number l2Error=0;
	number maxErr=0;

	bool bRes = true;
	// calculate scv size
	if (! calculate_vertex_vol(numsol,aaVolume)) {UG_LOG("ERROR: gradient computation failed in compute_error function!"); };

	//UG_LOG("----------------------------\n");

	if(!bRes) {UG_LOG("Error while calculating CV Volume.\n"); return false;}
	for (int si=0;si<numsol.num_subsets();++si)
	{
		// UG_LOG("*** " << si << "\n");
		for(VertexConstIterator iter = numsol.template begin<Vertex>(si);
									   iter != numsol.template end<Vertex>(si); ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
			number exactVal;
			position_accessor_type aaPos = domain.position_accessor();

		//	get vector holding all indices on the vertex
			std::vector<DoFIndex> ind;

			const size_t numInd = numsol.dof_indices(vrt, 0, ind);

		//	check indices
			if(numInd != 1) {UG_LOG("ERROR: Wrong number of indices!"); return false;}

			(*m_imDirichlet)(&exactVal,&aaPos[vrt],m_time,si,1);
			number differ = std::abs(DoFRef(numsol, ind[0])-exactVal);
		
			l1Error += aaVolume[vrt] * differ;
			l2Error += aaVolume[vrt] * differ*differ;

			if (m_print)
			{
			 //   if (differ>0)
			    UG_LOG("coord=" << aaPos[vrt] << " value=" << DoFRef(numsol, ind[0]) << " exact=" << exactVal << " error=" << differ << "\n");
			};
		
			if (differ > maxErr) maxErr = differ;

	     }
	};
	l2Error = sqrt(l2Error);
	UG_LOG("timestep " << m_timestep_nr << " time " << m_time << "\n");
	UG_LOG("l1 error: " << l1Error << "\n");
	UG_LOG("l2 error: " << l2Error << "\n");
	UG_LOG("maximum error: " << maxErr << "\n");
	return true;	
};

// next functions used for extrapolation equations as described in
// T.D. Aslam - A partial differential equation approach to multidimensional extrapolation JCP 193 2003

// compute gradient in vertices and volume of control volume in region given by sign of level set function
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::
calculate_vertex_grad_vol_sign(TGridFunction& u, aaGrad& aaGradient,aaVol& aaVolume,TGridFunction& phi,int sign)
{
	//	get domain
	domain_type& domain = *u.domain().get();

	//	get grid of domain
	typename domain_type::grid_type& grid = *domain.grid();

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	//	create a FV Geometry for the dimension
	DimFV1Geometry<dim> geo;

	//	hard code function (fct=0)
	//\todo: generalize
	size_t fct=0;

	// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	sum up all contributions of the sub control volumes to one vertex in an attachment
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)) continue;
		if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)) continue;
		if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)) continue;
		//	get iterators
		ElemIterator iter = u.template begin<ElemType>(si);
		ElemIterator iterEnd = u.template end<ElemType>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;

			//	get position accessor
			const position_accessor_type& aaPos = domain.position_accessor();

			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i)
			{
				vVrt[i] = elem->vertex(i);
				coCoord[i] = aaPos[vVrt[i]];
			};

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			//UG_LOG("Num Verts loaded: "<<vVrt.size()<<"\n");
			//UG_LOG("Num SCV computed: "<<geo.num_scv()<<"\n");

			number uValue[domain_traits<dim>::MaxNumVerticesOfElem];

			size_t noc = geo.num_scv();
			bool rightsign=true;
			for (size_t i=0;i < noc;i++)
			{
				//	get indices of function fct on vertex
				u.inner_dof_indices(vVrt[i], fct, multInd);

				//	read value of index from vector
				uValue[i]=DoFRef(u, multInd[0]);
				if (sign==-1)
				{
					if (DoFRef(phi, multInd[0])>0)
					{
						rightsign = false;
						break;
					}
				};
				if (sign==1)
				{
					if (DoFRef(phi, multInd[0])>0)
					{
						rightsign = false;
						break;
					};
				};
				//	debug log
				//UG_LOG("corner " << i << " " << uValue[i] << "\n");
			}
			if (! rightsign) continue;
			for (size_t i=0;i < noc;i++)
			{
				//	get indices of function fct on vertex
				u.inner_dof_indices(vVrt[i], fct, multInd);

				//	read value of index from vector
				uValue[i]=DoFRef(u, multInd[0]);

				//	debug log
				//UG_LOG("corner " << i << " " << uValue[i] << "\n");
			}

			//	storage for global gradient
			MathVector<dim> globalGrad;

			//	loop corners
			for (size_t i=0;i < noc;i++)
			{
				//	get scv for sh
				const typename DimFV1Geometry<dim>::SCV& scv = geo.scv(i);

				//	debug log
				//UG_LOG("gradient for corner " << i << "\n");

				//	reset global gradient
				globalGrad = 0.0;

				//	sum up gradients of shape functions in corner
				for(size_t sh = 0 ; sh < noc; ++sh)
				{
					//UG_LOG("local grad " << sh << " : " << scv.local_grad(sh) << "\n");
					//UG_LOG("unscaled global grad " << sh << " = " << scv.global_grad(sh) << "\n");
					//UG_LOG("uvalue(" << sh << ") =" << uValue[sh] << "\n");

					VecScaleAppend(globalGrad, uValue[sh], scv.global_grad(sh));
				}

				//	volume of scv
				number vol = scv.volume();

				//UG_LOG("*** global grad " << i << ": " << globalGrad << "\n");

				//	scale gradient by volume
				globalGrad *= vol;

				//	add both values to attachements
				aaGradient[vVrt[i]] += globalGrad;
				aaVolume[vVrt[i]] += vol;
			};
		}
	}

	// int count=0;
	position_accessor_type aaPos = u.domain()->position_accessor();
	for (int si=0;si < u.num_subsets();++si)
	{
		//UG_LOG("si " << si << "\n");
		for(VertexConstIterator iter = u.template begin<Vertex>(si);
			iter != u.template end<Vertex>(si); ++iter)
		{
			//	get vertex
			Vertex* vrt = *iter;
			if (aaVolume[vrt]!=0){
				(aaGradient[vrt]) /= aaVolume[vrt];
			};
			//exact[0] = cos(coord[0]);// 6*coord[0];
			//exact[1] = -4*sin(coord[1]);//-4*coord[1];
			//number gError = sqrt( (exact[0]-aaGradient[vrt][0])*(exact[0]-aaGradient[vrt][0]) + (exact[1]-aaGradient[vrt][1])*(exact[1]-aaGradient[vrt][1]) );
			//UG_LOG(count << "[ " << coord[0] << "," << coord[1] << " ] vol= " << aaVolume[vrt] << " " << "grad= ["
			//		<< aaGradient[vrt][0] << "," << aaGradient[vrt][1] << "] exact grad = ["
			//		<< exact[0] << "," << exact[1] << "] error: " << gError <<  "\n");
			// count++;
		}
	}
	
	//UG_LOG("#*#*#*#\n");
	return true;
}

// compute normal given by \frac{\nabla \phi}{|\nabla \phi|}
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_normal(TGridFunction& vx,TGridFunction& vy,TGridFunction& u)
{
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();

	//	create Attachment for scv-volume size
	ANumber aScvVolume;

	//	typedef of gradient attachment
	typedef Attachment<MathVector<dim> > AGradient;

	//	create Attachment for gradient
	AGradient aGradient;

	//	attach to grid
	grid.attach_to_vertices(aScvVolume);
	grid.attach_to_vertices(aGradient);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaVolume(grid, aScvVolume);
	Grid::VertexAttachmentAccessor<AGradient> aaGradient(grid, aGradient);

	// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

    MathVector<dim> coord;
    position_accessor_type aaPos = u.domain()->position_accessor();

    //	get vector holding all indices on the vertex
	std::vector<DoFIndex> ind;

    //	read indices on vertex
    // calculate scv size and gradient
    if (!calculate_vertex_grad_vol(u,aaGradient, aaVolume)) {UG_LOG("ERROR: gradient computation failed!"); return false;};
	if (m_limiter)
		limit_grad(u,aaGradient);
    // SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0); for debug set gradient to 0
	for (int si=0;si<u.num_subsets();++si)
	{
	    if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)) continue;
	    if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)) continue;
	    if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)) continue;
        VertexConstIterator iter = u.template begin<Vertex>(si);
		VertexConstIterator iterEnd = u.template end<Vertex>(si);
		for (;iter != iterEnd; ++iter)
		{
		    Vertex* vrt = *iter;
			u.inner_dof_indices(vrt, 0, ind);
			number vnorm = VecLength(aaGradient[vrt]);
			if (vnorm>1e-15)
			{
	    		DoFRef(vx, ind[0]) = aaGradient[vrt][0]/vnorm;
	    		DoFRef(vy, ind[0]) = aaGradient[vrt][1]/vnorm;
			} else {
				DoFRef(vx, ind[0]) = 0;
				DoFRef(vy, ind[0]) = 0;
			};
		};
    };
	return true;
}

// compute directional derivative in normal direction given by \normal \cdot \nabla u (see Aslam p. 2)
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_dnormal(TGridFunction& dnormal,TGridFunction& vx,TGridFunction& vy,TGridFunction& phi,TGridFunction& u)
{
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();

	//	create Attachment for scv-volume size
	ANumber aScvVolume;

	//	typedef of gradient attachment
	typedef Attachment<MathVector<dim> > AGradient;

	//	create Attachment for gradient
	AGradient aGradient;

	//	attach to grid
	grid.attach_to_vertices(aScvVolume);
	grid.attach_to_vertices(aGradient);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaVolume(grid, aScvVolume);
	Grid::VertexAttachmentAccessor<AGradient> aaGradient(grid, aGradient);

	// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

	MathVector<dim> coord;

	//	get vector holding all indices on the vertex
	std::vector<DoFIndex> ind;
	//	read indices on vertex
	position_accessor_type aaPos = u.domain()->position_accessor();
//UG_LOG("-------------\n");
	// calculate scv size and gradient of u
    if (! calculate_vertex_grad_vol_sign(u,aaGradient, aaVolume,phi,-1)) {UG_LOG("ERROR: gradient computation failed!"); return false;};
//	if (m_limiter)
//		limit_grad(u,aaGradient);
	// if (calculate_vertex_grad_vol(u,aaGradient, aaVolume)==false){UG_LOG("ERROR: gradient computation failed!"); return false;};
	// calculate normal of phi
	compute_normal(vx,vy,phi);
	for (int si=0;si<u.num_subsets();++si)
	{
		VertexConstIterator iter = u.template begin<Vertex>(si);
	    VertexConstIterator iterEnd = u.template end<Vertex>(si);
//	    UG_LOG("START INDEX" << si << "\n");
		for (;iter != iterEnd; ++iter)
		{
			Vertex* vrt = *iter;
			u.inner_dof_indices(vrt, 0, ind);
			coord = aaPos[vrt];
			DoFRef(dnormal, ind[0]) = DoFRef(vx, ind[0]) * aaGradient[vrt][0] + DoFRef(vy, ind[0]) * aaGradient[vrt][1];
//			if (DoFRef(phi, ind[0]) < 0)
//			    UG_LOG("coord=(" << coord[0] << "," << coord[1] << ") exact=" << coord[0]/sqrt(coord[0]*coord[0]+coord[1]*coord[1]) << " dnormal=" << BlockRef(dnormal[ind[0][0]],ind[0][1]) << "\n" << " v=(" << BlockRef(vx[ind[0][0]],ind[0][1]) << "," << BlockRef(vy[ind[0][0]],ind[0][1]) << ") " << " grad=(" << aaGradient[vrt][0] << "," << aaGradient[vrt][1] << ")" << "\n");
	    }
//		UG_LOG("#\n");
	};
	return true;
};

// compute directional derivative in normal direction of directional derivative in normal direction given by \normal \cdot \nabla (\normal \cdot \nabla u) (see Aslam p. 3)
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_ddnormal(TGridFunction& ddnormal,TGridFunction& dnormal,TGridFunction& vx,TGridFunction& vy,TGridFunction& phi,TGridFunction& u)
{
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();

	//	create Attachment for scv-volume size
	ANumber aScvVolume;

	//	typedef of gradient attachment
	typedef Attachment<MathVector<dim> > AGradient;

	//	create Attachment for gradient
	AGradient aGradient;

	//	attach to grid
	grid.attach_to_vertices(aScvVolume);
	grid.attach_to_vertices(aGradient);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaVolume(grid, aScvVolume);
	Grid::VertexAttachmentAccessor<AGradient> aaGradient(grid, aGradient);

	// initialize attachment value
	SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
	SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

	MathVector<dim> coord;

	//	read indices on vertex
	position_accessor_type aaPos = u.domain()->position_accessor();

	// calculate scv size and gradient of u
	compute_dnormal(dnormal,vx,vy,phi,u);
    if (! calculate_vertex_grad_vol_sign(dnormal,aaGradient, aaVolume,phi,-1)) {UG_LOG("ERROR: gradient computation failed!"); return false;};
    //if (m_limiter)
    //	limit_grad(dnormal,aaGradient);
	// if (! calculate_vertex_grad_vol(u,aaGradient, aaVolume)) {UG_LOG("ERROR: gradient computation failed!"); return false;};
	// calculate normal of phi

	std::vector<DoFIndex> ind;

	for (int si=0;si<u.num_subsets();++si)
	{
		VertexConstIterator iter = u.template begin<Vertex>(si);
	    VertexConstIterator iterEnd = u.template end<Vertex>(si);
////	    UG_LOG("START INDEX" << si << "\n");
		for (;iter != iterEnd; ++iter)
		{
			Vertex* vrt = *iter;
			u.inner_dof_indices(vrt, 0, ind);
			coord = aaPos[vrt];
			DoFRef(ddnormal, ind[0])
				= DoFRef(vx, ind[0]) * aaGradient[vrt][0]
				+ DoFRef(vy, ind[0]) * aaGradient[vrt][1];
		//	if (DoFRef(phi, ind[0])<0)
		//	    UG_LOG("coord=(" << coord[0] << "," << coord[1] << ") exact=" << coord[0]/sqrt(coord[0]*coord[0]+coord[1]*coord[1]) << " ddnormal=" << BlockRef(dnormal[ind[0][0]],ind[0][1]) << "\n" << " v=(" << BlockRef(vx[ind[0][0]],ind[0][1]) << "," << BlockRef(vy[ind[0][0]],ind[0][1]) << ") " << " grad=(" << aaGradient[vrt][0] << "," << aaGradient[vrt][1] << ")" << "\n");
		}
///		UG_LOG("#\n");
	};
	return true;
};

// in region given by level set sign:
// overwrite values in parameter unew with values of uold
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::overwrite
(
	TGridFunction& unew,
	TGridFunction& uold,
	TGridFunction& phi,
	int sign
)
{
	for (int si=0;si<unew.num_subsets();++si)
	{
		for(VertexConstIterator iter = unew.template begin<Vertex>(si);
									   iter != unew.template end<Vertex>(si); ++iter)
		{
			Vertex* vrt = * iter;

		//	read indices on vertex
			std::vector<DoFIndex> ind;
			unew.inner_dof_indices(vrt, 0, ind);
		    number phiValue = DoFRef(phi, ind[0]);
			int nodeSign;
			if (phiValue<0)  nodeSign =-1;
			if (phiValue>0)  nodeSign = 1;
			if (phiValue==0) nodeSign = 0;
		    if (nodeSign==sign)	DoFRef(unew, ind[0]) = DoFRef(uold, ind[0]);
		}
	};
    return true;
};

// in region given by level set sign:
// overwrite values in parameter unew with parameter value
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::overwrite
(
	TGridFunction& unew,
	number value,
	TGridFunction& phi,
	int sign
)
{
	for (int si=0;si<unew.num_subsets();++si)
	{
		for(VertexConstIterator iter = unew.template begin<Vertex>(si);
									   iter != unew.template end<Vertex>(si); ++iter)
		{
			Vertex* vrt = * iter;

		//	read indices on vertex
			std::vector<MultiIndex<2> > ind;
			unew.inner_dof_indices(vrt, 0, ind);
		    number phiValue = DoFRef(phi, ind[0]);
			int nodeSign;
			if (phiValue<0)  nodeSign =-1;
			if (phiValue>0)  nodeSign = 1;
			if (phiValue==0) nodeSign = 0;
		    if (nodeSign==sign)	DoFRef(unew, ind[0]) = value;
		}
	};
    return true;
};

// assign subsets as given by level set function 
// ug subset system would have to be changed for this to be useful
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::update_ls_subsets(TGridFunction& phi)
{
	//	get domain of grid function
    domain_type& domain = *phi.domain().get();

    //	create Multiindex
    std::vector<DoFIndex> ind;

	//	get element iterator type
	m_inactive_sg.set_subset_handler(domain.subset_handler());
    for (int si=0;si<domain.subset_handler()->num_subsets();++si)
    {
    	UG_LOG("******************* si " << si << " **********************\n");
		if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)) continue;
        if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)) continue;
		ElemIterator iter = phi.template begin<ElemType>(si);
     	ElemIterator iterEnd = phi.template end<ElemType>(si);
     	//	loop elements of dimension
     	for(  ;iter !=iterEnd; ++iter)
     	{
  	    //	get Elem
        	ElemType* elem = *iter;

        	//	coord and vertex array
        	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

        	//	get vertices and extract corner coordinates
        	const size_t numVertices = elem->num_vertices();
        	for(size_t i = 0; i < numVertices; ++i)
        	{
        		vVrt[i] = elem->vertex(i);
        	};

	        //	resize corners
	        std::vector<MathVector<dim> > coCoord;
     		//	compute center of mass
     		MathVector<dim> center;
     		std::vector<MathVector<dim> > grad;
     		center=0;
      		int noc=elem->num_vertices();
      		std::vector<number> phiCo(noc);
           	for(int i = 0; i < noc; ++i)
           	{
	        	phi.inner_dof_indices(vVrt[i], 0, ind);
     			phiCo[i]=DoFRef(phi, ind[0]);
      		};
			int firstNonzero=-1;
			for (int j=0;j<noc;j++)
			{
			    if (phiCo[j]!=0)
			    {
				    firstNonzero=j;
				    break;
				};
			};
			if (firstNonzero==-1) // all element nodes are on zero ls
			     domain.subset_handler()->assign_subset(elem,m_inside_elements_si);
			// add nodes on zero ls to zero ls node subsetgroup
			for (int i=0;i<noc;i++)
			{
				if (phiCo[i]==0)
				{
    			    int oldindex = domain.subset_handler()->get_subset_index(vVrt[i]);
	            	if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(oldindex)) continue;
		            if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(oldindex)) continue;
                    domain.subset_handler()->assign_subset(vVrt[i],m_onls_nodes_si);
 				};
			};
			bool onls=false;
			for (int i=firstNonzero+1;i<noc;i++)
			{
			    if (phiCo[firstNonzero]*phiCo[i]<0)
			    {
			    	onls=true;
			    	break;
			    }
			};
			if (! onls)
			{
			    if (phiCo[firstNonzero]<0)
			    {
			    //	UG_LOG("si before " << domain.subset_handler().get_subset_index(elem));
				    domain.subset_handler()->assign_subset(elem,m_inside_elements_si);
				    UG_LOG("element is inside \n");
				//    UG_LOG("si after " << domain.subset_handler().get_subset_index(elem) << "\n");
				//    UG_LOG("-- nr of subsets: " << phi.num_subsets() << " " << m_inside_elements_si <<  "\n");
				};
				if (phiCo[firstNonzero]>0)
				{
				//	UG_LOG("si before " << domain.subset_handler().get_subset_index(elem));
				    domain.subset_handler()->assign_subset(elem,m_outside_elements_si);
				//    UG_LOG("si after " << domain.subset_handler().get_subset_index(elem) <<  "\n");
				//    UG_LOG("|| nr of subsets: " << phi.num_subsets() << " " << m_outside_elements_si << "\n");
				};
			};
			if (onls)
			{
				//UG_LOG("si before " << domain.subset_handler().get_subset_index(elem));
			    domain.subset_handler()->assign_subset(elem,m_onls_elements_si);
			    //UG_LOG("si after " << domain.subset_handler().get_subset_index(elem) <<  "\n");
			    //UG_LOG("// nr of subsets: " << phi.num_subsets() << " " << m_onls_elements_si << "\n");
			    for (int i=0;i<noc;i++)
			    {
   					int oldindex = domain.subset_handler()->get_subset_index(vVrt[i]);
	    		    if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(oldindex)) continue;
       	            if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(oldindex)) continue;
				    if (phiCo[i]<0)
				    {
					    domain.subset_handler()->assign_subset(vVrt[i],m_inside_nodes_si);
					    UG_LOG("node is inside \n");
					};
					if (phiCo[i]>0)
					{
					    domain.subset_handler()->assign_subset(vVrt[i],m_outside_nodes_si);
					};
					if (phiCo[i]==0)
					{
						domain.subset_handler()->assign_subset(vVrt[i],m_onls_nodes_si);
					}
				};
			};
	    };
	};
/*  for(int sindex = 0; sindex < phi.num_subsets(); ++sindex)
    {
       	 UG_LOG("si " << sindex << "\n");
    	  //	get iterators
    	 ElemIterator iter = phi.template begin<ElemType>(sindex);
    	 ElemIterator iterEnd = phi.template end<ElemType>(sindex);
    	 int count=0;
    	for(  ;iter !=iterEnd; ++iter)
    	 {
    		ElemType* elem = *iter;
    		++count;
    	};
    	UG_LOG(count << " elements in subset\n");
    	count=0;
       	 for(VertexConstIterator iter = phi.template begin<Vertex>(sindex);
    	   iter != phi.template end<Vertex>(sindex); ++iter)
       	 {
       	     Vertex* vrt = *iter;
       	    ++count;
       	 }
    	 UG_LOG(count << " nodes in subset\n");
    };*/
  	return true;
}

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::init_ls_subsets(TGridFunction& phi)
{
    create_ls_subsets(phi);
	if (! update_ls_subsets(phi)) return false;
    return true;
}

template<typename TGridFunction>
void FV1LevelSetDisc<TGridFunction>::create_ls_subsets(TGridFunction& phi)
{
	//	get domain
	domain_type& domain = *phi.domain().get();
	UG_LOG("nr of subsets: " << domain.subset_handler()->num_subsets() << "\n");
   	m_inside_elements_si = domain.subset_handler()->num_subsets();
   	m_outside_elements_si = m_inside_elements_si + 1;
   	m_onls_elements_si = m_inside_elements_si + 2;
   	m_inside_nodes_si = m_inside_elements_si + 3;
   	m_outside_nodes_si = m_inside_elements_si + 4;
   	m_onls_nodes_si = m_inside_elements_si + 5;
   	domain.subset_handler()->subset_required(phi.num_subsets()+5);
   	UG_LOG("nr of subsets:" << domain.subset_handler()->num_subsets() << "\n");
};

// for runtime testing, delete later
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::runtimetest(TGridFunction& uNew)
{
	//	get domain of grid function
	domain_type& domain = *uNew.domain().get();

 	//	create a FV Geometry for the dimension
	DimFV1Geometry<dim> geo;

//	get position accessor
	const position_accessor_type& aaPos = domain.position_accessor();

//	resize corners
//	std::vector<MathVector<dim> > coCoord;
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	for (int si=0;si<uNew.num_subsets();++si)
	{
	    //	get iterators
	    ElemIterator iter = uNew.template begin<ElemType>(si);
	    ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
		//	get Elem
			ElemType* elem = *iter;

		//	extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i)
				coCoord[i] = aaPos[elem->vertex(i)];

		//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
		};
	};

	/*
	for (int si=0;si<uNew.num_subsets();++si)
	{
		//	get iterators
		ElemIterator iter = uNew.template begin<ElemType>(si);
		ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;
			//	get vertices of the Elem
			std::vector<Vertex*> vVrt;
			CollectVertices(vVrt, grid, elem);
			//	get position accessor
			const position_accessor_type& aaPos = domain.position_accessor();
			//	resize corners
			std::vector<MathVector<dim> > coCoord;
			//	extract corner coordinates
			for(size_t i = 0; i < vVrt.size(); ++i)
				coCoord.push_back( aaPos[vVrt[i]] );
		};
	};
	

	for (int si=0;si<domain.subset_handler().num_subsets();++si)
	{
		//	get iterators
		ElemIterator iter = uNew.template begin<ElemType>(si);
		ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;
			//	get vertices of the Elem
			std::vector<Vertex*> vVrt;
			CollectVertices(vVrt, grid, elem);
		};
	};

	MathVector<dim> coord;
	position_accessor_type aaPos = domain.position_accessor();

	for (int si=0;si<uNew.num_subsets();++si)
	{
		VertexConstIterator iter = uNew.template begin<Vertex>(si);
	    VertexConstIterator iterEnd = uNew.template end<Vertex>(si);
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			Vertex* vrt = *iter;
			coord = aaPos[vrt];
		}
	};*/
	return true;
}

} // end namespace LevelSet
} // end namespace ug

#endif /* LEVEL_SET_UTIL_IMPL_H_ */
