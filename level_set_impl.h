/*
 * level_set_util.h
 *
 *  Created on: 01.07.2011
 *      Author: Christian Wehner  christian.wehner@gcsc.uni-frankfurt.de
 */

#ifndef LEVEL_SET_UTIL_IMPL_H_
#define LEVEL_SET_UTIL_IMPL_H_

#include "level_set.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_grid/algorithms/attachment_util.h"

#include <algorithm>

namespace ug{
namespace LevelSet{

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::analytic_velocity(MathVector<dim> & v,number t,MathVector<dim> x)
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
number FV1LevelSetDisc<TGridFunction>::analytic_source(number t,MathVector<dim> x)
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
number FV1LevelSetDisc<TGridFunction>::analytic_solution(number t,MathVector<dim> x)
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
      //      z=min(min(abs(x[0]-1),abs(x[0]+1)),min(abs(x[1]-1),abs(x[1]+1)));
            //z=min(min(abs(x[0]),abs(1-x[0])),min(abs(x[1]),abs(1-x[1])));
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
bool FV1LevelSetDisc<TGridFunction>::fill_v_vec(TGridFunction& vel,int component)
{
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

	if (component>dim){
		UG_THROW("fill_v_vec: component > dim.");
	}

	position_accessor_type aaPos = vel.domain()->position_accessor();
	for (int si=0;si<vel.num_subsets();++si)
	{
		for(VertexBaseConstIterator iter = vel.template begin<VertexBase>(si);
									   iter != vel.template end<VertexBase>(si); ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			MathVector<dim> coord;
			MathVector<dim> vnode;
			coord = aaPos[vrt];

		//	get vector holding all indices on the vertex
			std::vector<MultiIndex<2> > ind;

			vel.inner_multi_indices(vrt, 0, ind);
			analytic_velocity(vnode,m_time,coord);
			BlockRef(vel[ind[0][0]],ind[0][1]) = vnode[component];
	     }
	};
	return true;
};

// next functions used for extrapolation equations as described in
// T.D. Aslam - A partial differential equation approach to multidimensional extrapolation JCP 193 2003

// compute normal given by \frac{\nabla \phi}{|\nabla \phi|}
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_normal(TGridFunction& vx,TGridFunction& vy,TGridFunction& u)
{
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

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
	std::vector<MultiIndex<2> > ind;

    //	read indices on vertex
    // calculate scv size and gradient
    if (calculate_vertex_grad_vol(u,aaGradient, aaVolume)==false){UG_LOG("ERROR: gradient computation failed!"); return false;};
	if (m_limiter==true){
		limit_grad(u,aaGradient);
	};
    // SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0); for debug set gradient to 0
	for (int si=0;si<u.num_subsets();++si){
	    if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
	    if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)==true) continue;
	    if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)==true) continue;
        VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
		VertexBaseConstIterator iterEnd = u.template end<VertexBase>(si);
		for (;iter != iterEnd; ++iter){
		    VertexBase* vrt = *iter;
			u.inner_multi_indices(vrt, 0, ind);
			number vnorm = VecLength(aaGradient[vrt]);
			if (vnorm>1e-15){
	    		BlockRef(vx[ind[0][0]],ind[0][1]) = aaGradient[vrt][0]/vnorm;
		    	BlockRef(vy[ind[0][0]],ind[0][1]) = aaGradient[vrt][1]/vnorm;
			} else {
			    BlockRef(vx[ind[0][0]],ind[0][1]) = 0;
    			BlockRef(vy[ind[0][0]],ind[0][1]) = 0;
			};
		};
    };
	return true;
}

// compute directional derivative in normal direction given by \normal \cdot \nabla u (see Aslam p. 2)
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_dnormal(TGridFunction& dnormal,TGridFunction& vx,TGridFunction& vy,TGridFunction& phi,TGridFunction& u)
{
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

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
	std::vector<MultiIndex<2> > ind;
	//	read indices on vertex
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type aaPos = u.domain()->position_accessor();
///UG_LOG("-------------\n");
	// calculate scv size and gradient of u
    if (calculate_vertex_grad_vol_sign(u,aaGradient, aaVolume,phi,-1)==false){UG_LOG("ERROR: gradient computation failed!"); return false;};
	if (m_limiter==true){
//		limit_grad(u,aaGradient);
	};
	// if (calculate_vertex_grad_vol(u,aaGradient, aaVolume)==false){UG_LOG("ERROR: gradient computation failed!"); return false;};
	// calculate normal of phi
	compute_normal(vx,vy,phi);
	for (int si=0;si<u.num_subsets();++si){
		VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
	    VertexBaseConstIterator iterEnd = u.template end<VertexBase>(si);
////	    UG_LOG("START INDEX" << si << "\n");
		for (;iter != iterEnd; ++iter){
			VertexBase* vrt = *iter;
			u.inner_multi_indices(vrt, 0, ind);
			coord = aaPos[vrt];
			BlockRef(dnormal[ind[0][0]],ind[0][1]) = BlockRef(vx[ind[0][0]],ind[0][1]) * aaGradient[vrt][0] + BlockRef(vy[ind[0][0]],ind[0][1]) * aaGradient[vrt][1];
///			if (BlockRef(phi[ind[0][0]],ind[0][1])<0)
///			    UG_LOG("coord=(" << coord[0] << "," << coord[1] << ") exact=" << coord[0]/sqrt(coord[0]*coord[0]+coord[1]*coord[1]) << " dnormal=" << BlockRef(dnormal[ind[0][0]],ind[0][1]) << "\n" << " v=(" << BlockRef(vx[ind[0][0]],ind[0][1]) << "," << BlockRef(vy[ind[0][0]],ind[0][1]) << ") " << " grad=(" << aaGradient[vrt][0] << "," << aaGradient[vrt][1] << ")" << "\n");
	    }
///		UG_LOG("#\n");
	};
	return true;
};

// compute directional derivative in normal direction of directional derivative in normal direction given by \normal \cdot \nabla (\normal \cdot \nabla u) (see Aslam p. 3)
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_ddnormal(TGridFunction& ddnormal,TGridFunction& dnormal,TGridFunction& vx,TGridFunction& vy,TGridFunction& phi,TGridFunction& u)
{
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

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
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type aaPos = u.domain()->position_accessor();

	// calculate scv size and gradient of u
	compute_dnormal(dnormal,vx,vy,phi,u);
    if (calculate_vertex_grad_vol_sign(dnormal,aaGradient, aaVolume,phi,-1)==false){UG_LOG("ERROR: gradient computation failed!"); return false;};
    if (m_limiter==true){
    //	limit_grad(dnormal,aaGradient);
    };
	// if (calculate_vertex_grad_vol(u,aaGradient, aaVolume)==false){UG_LOG("ERROR: gradient computation failed!"); return false;};
	// calculate normal of phi

	std::vector<MultiIndex<2> > ind;

	for (int si=0;si<u.num_subsets();++si){
		VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
	    VertexBaseConstIterator iterEnd = u.template end<VertexBase>(si);
////	    UG_LOG("START INDEX" << si << "\n");
		for (;iter != iterEnd; ++iter){
			VertexBase* vrt = *iter;
			u.inner_multi_indices(vrt, 0, ind);
			coord = aaPos[vrt];
			BlockRef(ddnormal[ind[0][0]],ind[0][1])
				= BlockRef(vx[ind[0][0]],ind[0][1]) * aaGradient[vrt][0]
				+ BlockRef(vy[ind[0][0]],ind[0][1]) * aaGradient[vrt][1];
		//	if (BlockRef(phi[ind[0][0]],ind[0][1])<0)
		//	    UG_LOG("coord=(" << coord[0] << "," << coord[1] << ") exact=" << coord[0]/sqrt(coord[0]*coord[0]+coord[1]*coord[1]) << " ddnormal=" << BlockRef(dnormal[ind[0][0]],ind[0][1]) << "\n" << " v=(" << BlockRef(vx[ind[0][0]],ind[0][1]) << "," << BlockRef(vy[ind[0][0]],ind[0][1]) << ") " << " grad=(" << aaGradient[vrt][0] << "," << aaGradient[vrt][1] << ")" << "\n");
		}
///		UG_LOG("#\n");
	};
	return true;
};

int bubblesort(std::vector<number> array,std::vector<int>& list){
    int i,j,mini,k;
    int length=array.size();
	number min;
	for (i=0;i<length;i++) list[i]=i;
	for (i=0;i<length-1;i++){
		min = array[i];
		mini= i;
		for (j=i+1;j<length;j++){
			if (array[j]<min){
				mini = j;
				min = array[j];
			};
		};
		array[mini] = array[i];
		k=list[i];
		list[i]     = list[mini];
		list[mini]  = k;
		array[i]    = min;
	};
	return 0;
};

/* Next functions are self-explaining linear algebra functions
   used below in function solveLS
*/
bool multMatVec(const std::vector<number>& avec, const std::vector<number>& b,std::vector<number>& c,size_t m,size_t n){
   size_t i;
   size_t j;
   size_t count=0;
   for(i = 0; i < m; i = i + 1){
	   c[i]=0.0;
       for(j = 0; j < n; j = j + 1){
            c[i] = c[i] + avec[count] * b[j];
		    count++;
	   }
   }
   return true;
}

/*
* Solve linear system for n x n - matrix given as array
* Remark: historic document because it was computed as exercise for
* "Algorithmische Optimierung I", WS 04/05
*/
bool solveLS(std::vector<number>& x,/* solution */
		const std::vector<number>& matField,/* matrix given as field */
            const std::vector<number>& b /* rhs */){
	size_t n=b.size();
	number P[n][n];
    std::vector<number> Pvec(n*n);
    std::vector<number> b2(n);
    number A[n][n];
	number L[n][n];
	number U[n][n];
	number z[n];
	size_t i,k,j,l;
	bool boolean=false;
	number d,max;
    size_t count=0;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            A[i][j]=matField[count];
            if (i!=j)
               P[i][j]=0;
            else
               P[i][j]=1;
            L[i][j]=0;
            U[i][j]=0;
            count ++;
        }
    }
	for (i=0;i<n;i++){
		j=i;
		// Suche groesstes Element in Spalten
		max=abs(A[i][i]);
		for (k=i+1;k<n;k++){
			if (abs(A[k][i])>max){
				j=k;
				max=abs(A[k][i]);
			};
		};
		//debug
		//UG_LOG(max << "\n");
		if (max<1e-14){
			// A nicht invertierbar
			// printf(".\n");
		    return false;
		};
		if (i!=j){
			// Vertauschen der Zeilen
			if (boolean==false){
				// beim ersten Durchlauf Vertauschen der Spalten der
				// Permutationsmatrix und der Zeilen von A,L.
				boolean=true;
				for (k=0;k<n;k++){
					d=P[k][i];
				    P[k][i]=P[k][j];
				    P[k][j]=d;
					d=A[i][k];
				    A[i][k]=A[j][k];
				    A[j][k]=d;
					d=L[i][k];
				    L[i][k]=L[j][k];
				    L[j][k]=d;
				};
			} else{
				for (k=0;k<n;k++){
				    // Vertauschen der Zeilen von P,A,L
				    d=P[i][k];
				    P[i][k]=P[j][k];
				    P[j][k]=d;
					d=A[i][k];
				    A[i][k]=A[j][k];
				    A[j][k]=d;
					d=L[i][k];
				    L[i][k]=L[j][k];
				    L[j][k]=d;
				};
			};
		};
		// Elimination
		L[i][i]=1;
		for (k=i+1;k<n;k++){
			L[k][i]=A[k][i]/A[i][i];
			for (l=i+1;l<n;l++){
				A[k][l]=A[k][l]-L[k][i]*A[i][l];
			};
		};
	};
	// U ist der obere Dreiecksteil von A
	for (i=0;i<n;i++){
		for (j=i;j<n;j++){
			U[i][j]=A[i][j];
		};
	};
	// A*x = b <-> P*A*x = P*b
    count = 0;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            Pvec[count]=P[i][j];
			count++;
        }
    }
	multMatVec(Pvec,b,b2,n,n);
	// Loese zuerst L*z = P*b
	for (i=0;i<n;i++){
		number s=b2[i];
		for (j=0;j<i;j++){
			s=s-z[j]*L[i][j];
		};
		z[i]=s;
	};
	// Loese U*x = z
	for (i=n-1;i>=0;i--){
		number s=z[i];
		for (j=n-1;j>i;j--){
			s=s-x[j]*U[i][j];
		};
		x[i]=s/U[i][i];
		if (i==0) break;
	};
	return true;
};


template<size_t n>
bool solveLS(MathVector<n>& x,const MathMatrix<n,n>& M,const MathVector<n>& b){
	std::vector<number> xvec(n);
	std::vector<number> bvec(n);
	std::vector<number> Mvec(n*n);
	for (size_t i=0;i<n;i++){
		bvec[i]=b[i];
		for (size_t j=0;j<n;j++){
			Mvec[n*i+j] = M[i][j];
		};
	}
	if (solveLS(xvec,Mvec,bvec)==false) return false;
	for (size_t i=0;i<n;i++) x[i] = xvec[i];
	return true;
}

template<size_t m,size_t n>
bool leastSquares(MathVector<n>& x,const MathMatrix<m,n>& M,const MathVector<m>& b){
	if(m<n){
		UG_THROW("Least squares method not suitable for m x n matrix with m<n.\n");
	}
	MathMatrix<n,n> M2;
	MathVector<n> b2;
	MatMultiplyMTM(M2,M);
	TransposedMatVecMult(b2,M,b);
	return solveLS(x,M2,b2);
}

bool leastSquares(std::vector<number>& x,const std::vector<number>& mField,const std::vector<number>& b){
	size_t m = b.size();
	size_t n = x.size();
	if (mField.size()!=m*n){
		n = mField.size()/b.size();
	}
	std::vector<number> tmmField(n*n);
	std::vector<number> tmb(n);
	number z;
	// compute A^t * A
	for (size_t i=0;i<n;i++){
		for (size_t j=i;j<n;j++){
			z=0;
			for (size_t k=0;k<m;k++){
				z+=mField[n*k+i]*mField[n*k+j];
			}
			tmmField[j*n+i]=z;
			if (j!=i) tmmField[i*n+j]=z;
		}
		tmb[i]=0;
		for (size_t k=0;k<m;k++) tmb[i]+= mField[n*k+i]*b[k];
	}
	return solveLS(x,tmmField,tmb);
}


/// averages positions by arithmetic mean
/**
 * Arithmetic Mean of Positions
 * returns the arithmetic mean of positions
 *
 * \param[in]  vCornerCoords	positions
 * \param[in]  num				number of positions
 * \param[out] vOut				arithmetic mean of positions

template <typename TPosition>
void FV1LevelSetDisc<TGridFunction>::AverageCoord(TPosition& vOut, const TPosition* vCornerCoords, size_t num)
{
	vOut = vCornerCoords[0];
	for(size_t j = 1; j < num; ++j)
	{
		vOut += vCornerCoords[j];
	}
	vOut *= 1./(number)num;
} */

/// averages positions by arithmetic mean
/**
 * Arithmetic Mean of Positions
 * returns the arithmetic mean of positions
 *
 * \param[in]  vCornerCoords	positions
 * \param[in]  num				number of positions
 * \param[out] vOut				arithmetic mean of positions
 */
template <typename TPosition>
void AverageCoords(TPosition& vOut, const std::vector<TPosition>& vCornerCoords, size_t num)
{
	vOut = vCornerCoords[0];
	for(size_t j = 1; j < num; ++j)
	{
		vOut += vCornerCoords[j];
	}
	vOut *= 1./(number)num;
}

// computation of curvature
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::computeElementCurvature2d(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order){
    std::vector<MathVector<dim> > interPoints;
    size_t nOfPoints=co.size();
    interPoints.resize(elementnoc);
    std::vector<number> distToBaseP(nOfPoints);
    std::vector<int> sortedList(nOfPoints);
    MathVector<dim> elemBasePoint;
    size_t interpointssize = 0;
    // compute element base point
    for (size_t i=0;i<elementnoc;i++){
        if (phi[i]==0){
            interPoints[interpointssize]=co[i];
            interpointssize++;
            continue;
        };
        for (size_t j=i+1;j<elementnoc;j++){
            if (phi[i]*phi[j]<0){
                interPoints[interpointssize]=co[i];
                interpointssize++;
            }
        }
    }
    AverageCoords(elemBasePoint,interPoints,interpointssize);
    //debugUG_LOG("base point: " << elemBasePoint << "\n");
    for (size_t i=0;i<nOfPoints;i++){
        distToBaseP[i]=VecDistance(co[i],elemBasePoint);
        //debugUG_LOG("dTBp[" << i << "]=" << distToBaseP[i] << "\n");
    };
    // sort nodes by distance to base point
    distToBaseP.resize(nOfPoints);
    sortedList.resize(nOfPoints);
    bubblesort(distToBaseP,sortedList);
    size_t criticalnofinterp;
    if (order==2) criticalnofinterp=5;
    if (order==3) criticalnofinterp=7;
    size_t nrOfInterPoints=0;
    std::vector<number> mat;
    size_t vlength=(size_t)round(0.5*(order+1)*(order+2));
    mat.resize(vlength*vlength);
	std::vector<number> interphi;
	std::vector<number> coeffs;
	interphi.resize(vlength);
	coeffs.resize(vlength);
	for (size_t i=criticalnofinterp*(vlength-1);i<vlength*vlength;i++){
        mat[i] = (number)(0.8*i*i-0.4*i)/(i+1)+i*i-0.5+sin(0.3*i)+cos(0.75*i)+i*i*i-exp(0.1*i);
    };
	size_t matindex=0;
	//debugUG_LOG("####### nofpoints=" << nOfPoints << "\n");
    for (size_t j=0;j<nOfPoints;j++){
        size_t i=sortedList[j];
        //debugUG_LOG("list[" << j << "]=" << sortedList[j] << "\n");
        for (size_t ii=0;ii<=order;ii++){
        	for (size_t k=0;k<=ii;k++){
        		//debugUG_LOG(i << " " << co[i] << "\n");
        		mat[matindex]=std::pow((number)co[i][0],(int)(ii-k))*std::pow((number)co[i][1],(int)k);
        		matindex++;
        	};
        };
        interphi[nrOfInterPoints]=phi[i];
        //debugUG_LOG("rhs[" << nrOfInterPoints << "] = " << interphi[nrOfInterPoints] << "\n");
        //debugfor (size_t kdebug=0;kdebug<nrOfInterPoints*vlength;kdebug++){
        	//debugUG_LOG(mat[kdebug] << "\n");
        //debug}
        //debugUG_LOG("------------------------------------------- " << j << "\n");
        if (j>=criticalnofinterp){
        	//debugUG_LOG("********************\n");
			if (solveLS(coeffs,mat,interphi)==false){
                matindex-=vlength;
                continue;
            };
        };
        //debugUG_LOG("co(" << nrOfInterPoints+1 << ",:)=[" << co[i][0] << "," << co[i][1] << "];" << " phi(" << nrOfInterPoints+1 << ")=" << phi[i] << ";\n");
        nrOfInterPoints++;
		if (nrOfInterPoints==vlength) break;
    };
	if (nrOfInterPoints<vlength){
		UG_LOG("warning: not enough interpolations points for desired order "<< order << " (found " << nrOfInterPoints << ", needed " << vlength << ") Reduce order.\n");
		if (order==1) return false;
		if (computeElementCurvature2d(kappa,elementnoc,co,phi,order-1)==true) return true; else return false;
	}
	number dxphi,dyphi,dxxphi,dxyphi,dyyphi;
	//debugUG_LOG("coeff.size =" << coeffs.size() << "\n");
	//debugfor (size_t j=0;j<coeffs.size();j++){
	//debugUG_LOG("coeff[" << j << "]=" << coeffs[j] << "\n");
	//debug}
	if (order==2){
        dxphi=coeffs[1]+2*coeffs[3]*elemBasePoint[0]+coeffs[4]*elemBasePoint[1];
        dyphi=coeffs[2]+coeffs[4]*elemBasePoint[0]+2*coeffs[5]*elemBasePoint[1];
    };
    if (order==3){
		dxphi=coeffs[1]+2*coeffs[3]*elemBasePoint[0]+coeffs[4]*elemBasePoint[1]+3*coeffs[6]*elemBasePoint[0]*elemBasePoint[0]+2*coeffs[7]*elemBasePoint[0]*elemBasePoint[1]+coeffs[8]*elemBasePoint[1]*elemBasePoint[1];
		dyphi=coeffs[2]+coeffs[4]*elemBasePoint[0]+2*coeffs[5]*elemBasePoint[1]+coeffs[7]*elemBasePoint[0]*elemBasePoint[0]+2*coeffs[8]*elemBasePoint[0]*elemBasePoint[1]+3*coeffs[9]*elemBasePoint[1]*elemBasePoint[1];
    };
    number gradnorm=sqrt(dxphi*dxphi+dyphi*dyphi);
	// characteristic element length
    number hh= std::min(VecDistance(co[0],co[1]),VecDistance(co[0],co[2]));
    number coelemBasePoint[dim];
    number sol[dim];
    coelemBasePoint[0] = elemBasePoint[0] + 0.5*hh*dxphi/gradnorm;
    coelemBasePoint[1] = elemBasePoint[1] + 0.5*hh*dyphi/gradnorm;
    /* find intersection of line between elemBasePoint and coelemBaseP with zero level set */
    number itgamma = 0;
    number phival,phigamma;
    size_t nrofiterations=0;
    for (;nrofiterations<50;nrofiterations++){
        number b = coelemBasePoint[0]-elemBasePoint[0];
        number d = coelemBasePoint[1]-elemBasePoint[1];
        sol[0] = elemBasePoint[0]+itgamma*(coelemBasePoint[0]-elemBasePoint[0]);
        sol[1] = elemBasePoint[1]+itgamma*(coelemBasePoint[1]-elemBasePoint[1]);
        if (order==2)
            phival = coeffs[0]+coeffs[1]*sol[0]+coeffs[2]*sol[1]+coeffs[3]*sol[0]*sol[0]+coeffs[4]*sol[0]*sol[1]+coeffs[5]*sol[1]*sol[1];
        if (order==3)
            phival = coeffs[0]+coeffs[1]*sol[0]+coeffs[2]*sol[1]+coeffs[3]*sol[0]*sol[0]+coeffs[4]*sol[0]*sol[1]+coeffs[5]*sol[1]*sol[1]
					+coeffs[6]*sol[0]*sol[0]*sol[0]+coeffs[7]*sol[0]*sol[0]*sol[1]+coeffs[8]*sol[0]*sol[1]*sol[1]+coeffs[9]*sol[1]*sol[1]*sol[1];
        if (abs(phival)<1e-12) break;
        phigamma = coeffs[1]*b+coeffs[2]*d+2*coeffs[3]*sol[0]*b+coeffs[4]*b*sol[1]+coeffs[4]*sol[0]*d+2*coeffs[5]*sol[1]*d+3*coeffs[6]*sol[0]*sol[0]*b
					+2*coeffs[7]*sol[0]*sol[1]*b+coeffs[7]*sol[0]*sol[0]*d+coeffs[8]*b*sol[1]*sol[1]+2*coeffs[8]*sol[0]*sol[1]*d+3*coeffs[9]*sol[1]*sol[1]*d;
        itgamma -= phival/phigamma;
    };
    if (nrofiterations==80){
        UG_THROW("Diverging Newton method in curvature computation (error=" << phival << ")\n");
        return false;
    };
    if (order==2){
        dxphi=coeffs[1]+2*coeffs[3]*sol[0]+coeffs[4]*sol[1];
        dyphi=coeffs[2]+coeffs[4]*sol[0]+2*coeffs[5]*sol[1];
        dxxphi=2*coeffs[3];
        dxyphi=coeffs[4];
        dyyphi=2*coeffs[5];
    };
    if (order==3){
		dxphi=coeffs[1]+2*coeffs[3]*sol[0]+coeffs[4]*sol[1]+3*coeffs[6]*sol[0]*sol[0]+2*coeffs[7]*sol[0]*sol[1]+coeffs[8]*sol[1]*sol[1];
		dyphi=coeffs[2]+coeffs[4]*sol[0]+2*coeffs[5]*sol[1]+coeffs[7]*sol[0]*sol[0]+2*coeffs[8]*sol[0]*sol[1]+3*coeffs[9]*sol[1]*sol[1];
        dxxphi=2*coeffs[3]+6*coeffs[6]*sol[0]+2*coeffs[7]*sol[1];
        dxyphi=coeffs[4]+2*coeffs[7]*sol[0]+2*coeffs[8]*sol[1];
        dyyphi=2*coeffs[5]+2*coeffs[8]*sol[0]+6*coeffs[9]*sol[1];
    };
    number t = sqrt(dxphi*dxphi+dyphi*dyphi);
    kappa = - (dyyphi * dxphi * dxphi - 2 * dxyphi * dxphi * dyphi + dxxphi * dyphi * dyphi) / (t * t * t);
//  UG_LOG("kappa=" << kappa << "\n");
//	kappa = -3.333333333333333333333333333333333333333;// for debug
//	kappa = -5;// for debug
    //debugUG_THROW("Grid Level Type not in ['top' | 'surf'].");
    return true;
}

// computation of curvature
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::computeElementCurvature2d2(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order,number nodefactor){
	if (order<2) return false;
    std::vector<MathVector<dim> > interPoints;
    size_t nOfPoints=co.size();
    interPoints.resize(elementnoc);
    std::vector<number> distToBaseP(nOfPoints);
    std::vector<int> sortedList(nOfPoints);
    MathVector<dim> elemBasePoint;
    size_t interpointssize = 0;
    // compute element base point
    for (size_t i=0;i<elementnoc;i++){
        if (phi[i]==0){
            interPoints[interpointssize]=co[i];
            interpointssize++;
            continue;
        };
        for (size_t j=i+1;j<elementnoc;j++){
            if (phi[i]*phi[j]<0){
                interPoints[interpointssize]=co[i];
                interpointssize++;
            }
        }
    }
    AverageCoords(elemBasePoint,interPoints,interpointssize);
    //debugUG_LOG("base point: " << elemBasePoint << "\n");
    for (size_t i=0;i<nOfPoints;i++){
        distToBaseP[i]=VecDistance(co[i],elemBasePoint);
        //debugUG_LOG("dTBp[" << i << "]=" << distToBaseP[i] << "\n");
    };
    // sort nodes by distance to base point
    distToBaseP.resize(nOfPoints);
    sortedList.resize(nOfPoints);
    bubblesort(distToBaseP,sortedList);
//	for (size_t i=0;i<nOfPoints;i++) sortedList[i]=i;
    size_t vlength=(size_t)round(0.5*(order+1)*(order+2));
	size_t nrOfInterPoints = round(vlength*nodefactor);
	size_t totalNrOfInterpoints = co.size();
	if (totalNrOfInterpoints<nrOfInterPoints){
		UG_LOG("Warning: not enough nodes for desired order " << order << " and least squares factor " << nodefactor << ". Needed " << nrOfInterPoints << " nodes, given " << totalNrOfInterpoints << " nodes. Reduce order to " << order-1 << ".\n");
		return computeElementCurvature2d2(kappa,elementnoc,co,phi,order-1,nodefactor);
	};
	std::vector<number> interM(vlength*nrOfInterPoints);
	std::vector<number> coeffs(vlength);
	std::vector<number> interRhs(nrOfInterPoints);
	size_t matindex = 0;
    for (size_t j=0;j<nrOfInterPoints;j++){
        size_t i=sortedList[j];
		// interpolation rhs
		interRhs[j] = phi[i];
		// interpolation matrix
        for (size_t ii=0;ii<=order;ii++){
        	for (size_t k=0;k<=ii;k++){
        		interM[matindex]=std::pow((number)co[i][0],(int)(ii-k))*std::pow((number)co[i][1],(int)k);
        		matindex++;
        	};
        };
    };
	if (leastSquares(coeffs,interM,interRhs)==false){
		UG_LOG("Least squares problem had no regular solution. Reduce order to " << order-1 << ".\n");
/*		for (size_t i=0;i<interM.size();i++){
			UG_LOG(interM[i] << " ");
		}
		UG_LOG("-----------\n");
		for (size_t i=0;i<nrOfInterPoints;i++){
			UG_LOG(interRhs[i] << " ");
		}
		UG_LOG("\n");*/
		return computeElementCurvature2d2(kappa,elementnoc,co,phi,order-1,nodefactor);
	};
	number dxphi,dyphi,dxxphi,dxyphi,dyyphi;
	if (order==2){
        dxphi=coeffs[1]+2*coeffs[3]*elemBasePoint[0]+coeffs[4]*elemBasePoint[1];
        dyphi=coeffs[2]+coeffs[4]*elemBasePoint[0]+2*coeffs[5]*elemBasePoint[1];
    };
    if (order==3){
		dxphi=coeffs[1]+2*coeffs[3]*elemBasePoint[0]+coeffs[4]*elemBasePoint[1]+3*coeffs[6]*elemBasePoint[0]*elemBasePoint[0]+2*coeffs[7]*elemBasePoint[0]*elemBasePoint[1]+coeffs[8]*elemBasePoint[1]*elemBasePoint[1];
		dyphi=coeffs[2]+coeffs[4]*elemBasePoint[0]+2*coeffs[5]*elemBasePoint[1]+coeffs[7]*elemBasePoint[0]*elemBasePoint[0]+2*coeffs[8]*elemBasePoint[0]*elemBasePoint[1]+3*coeffs[9]*elemBasePoint[1]*elemBasePoint[1];
    };
    number gradnorm=sqrt(dxphi*dxphi+dyphi*dyphi);
	// characteristic element length
    number hh= min(VecDistance(co[0],co[1]),VecDistance(co[0],co[2]));
    number coelemBasePoint[dim];
    number sol[dim];
    coelemBasePoint[0] = elemBasePoint[0] + 0.5*hh*dxphi/gradnorm;
    coelemBasePoint[1] = elemBasePoint[1] + 0.5*hh*dyphi/gradnorm;
    /* find intersection of line between elemBasePoint and coelemBaseP with zero level set */
    number itgamma = 0;
    number phival,phigamma;
    size_t nrofiterations=0;
    for (;nrofiterations<50;nrofiterations++){
        number b = coelemBasePoint[0]-elemBasePoint[0];
        number d = coelemBasePoint[1]-elemBasePoint[1];
        sol[0] = elemBasePoint[0]+itgamma*(coelemBasePoint[0]-elemBasePoint[0]);
        sol[1] = elemBasePoint[1]+itgamma*(coelemBasePoint[1]-elemBasePoint[1]);
        if (order==2)
            phival = coeffs[0]+coeffs[1]*sol[0]+coeffs[2]*sol[1]+coeffs[3]*sol[0]*sol[0]+coeffs[4]*sol[0]*sol[1]+coeffs[5]*sol[1]*sol[1];
        if (order==3)
            phival = coeffs[0]+coeffs[1]*sol[0]+coeffs[2]*sol[1]+coeffs[3]*sol[0]*sol[0]+coeffs[4]*sol[0]*sol[1]+coeffs[5]*sol[1]*sol[1]
					+coeffs[6]*sol[0]*sol[0]*sol[0]+coeffs[7]*sol[0]*sol[0]*sol[1]+coeffs[8]*sol[0]*sol[1]*sol[1]+coeffs[9]*sol[1]*sol[1]*sol[1];
        if (abs(phival)<1e-12) break;
        phigamma = coeffs[1]*b+coeffs[2]*d+2*coeffs[3]*sol[0]*b+coeffs[4]*b*sol[1]+coeffs[4]*sol[0]*d+2*coeffs[5]*sol[1]*d+3*coeffs[6]*sol[0]*sol[0]*b
					+2*coeffs[7]*sol[0]*sol[1]*b+coeffs[7]*sol[0]*sol[0]*d+coeffs[8]*b*sol[1]*sol[1]+2*coeffs[8]*sol[0]*sol[1]*d+3*coeffs[9]*sol[1]*sol[1]*d;
        itgamma -= phival/phigamma;
    };
    if (nrofiterations==80){
        UG_THROW("Diverging Newton method in curvature computation (error=" << phival << ")\n");
        return false;
    };
    if (order==2){
        dxphi=coeffs[1]+2*coeffs[3]*sol[0]+coeffs[4]*sol[1];
        dyphi=coeffs[2]+coeffs[4]*sol[0]+2*coeffs[5]*sol[1];
        dxxphi=2*coeffs[3];
        dxyphi=coeffs[4];
        dyyphi=2*coeffs[5];
    };
    if (order==3){
		dxphi=coeffs[1]+2*coeffs[3]*sol[0]+coeffs[4]*sol[1]+3*coeffs[6]*sol[0]*sol[0]+2*coeffs[7]*sol[0]*sol[1]+coeffs[8]*sol[1]*sol[1];
		dyphi=coeffs[2]+coeffs[4]*sol[0]+2*coeffs[5]*sol[1]+coeffs[7]*sol[0]*sol[0]+2*coeffs[8]*sol[0]*sol[1]+3*coeffs[9]*sol[1]*sol[1];
        dxxphi=2*coeffs[3]+6*coeffs[6]*sol[0]+2*coeffs[7]*sol[1];
        dxyphi=coeffs[4]+2*coeffs[7]*sol[0]+2*coeffs[8]*sol[1];
        dyyphi=2*coeffs[5]+2*coeffs[8]*sol[0]+6*coeffs[9]*sol[1];
    };
    number t = sqrt(dxphi*dxphi+dyphi*dyphi);
    kappa = - (dyyphi * dxphi * dxphi - 2 * dxyphi * dxphi * dyphi + dxxphi * dyphi * dyphi) / (t * t * t);
    return true;
}


template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::computeElementCurvatureOnGrid2d(TGridFunction& u,size_t order,number leastSquaresFactor){
	typedef typename TGridFunction::domain_type domain_type;
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& aaPos = u.domain()->position_accessor();
	std::vector<MultiIndex<2> > ind;
	number maxnormerr = 0;
	//	loop elements of dimension
	for (int si=0;si<u.num_subsets();++si){
		ElemIterator iter = u.template begin<ElemType>(si);
		ElemIterator iterEnd = u.template end<ElemType>(si);
		if (u.num_fct(si)<2){
			UG_THROW("No curvature component in approximation space.");
		}
		if (u.local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, 1)){
			UG_THROW("First component in approximation space must be of Lagrange 1 type.");
		}
		if (u.local_finite_element_id(1) != LFEID(LFEID::PIECEWISE_CONSTANT,0)){
			UG_THROW("Second component in approximation space must be of piecewise constant type.");
		}
		std::vector<MathVector<dim> > coord;
		std::vector<number> phi;
		std::vector<VertexBase*> nbrs;
		std::vector<VertexBase*> nbrCandidates;
		size_t depth=order;
		std::vector<size_t> stageStart(depth+2);
	for(  ;iter !=iterEnd; ++iter)
	{
	    //	get Elem
		ElemType* elem = *iter;
		//	get position accessor
		size_t noc=elem->num_vertices();
		nbrs.resize(noc);
		phi.resize(noc);
		// check if zero level set is on element
		bool onls=false;
		bool nonzerophifound=false;
		number nonzerophi;
		//debug		UG_LOG("noc = " << noc << "\n");
		for (size_t i=0;i<noc;i++){
			nbrs[i]=elem->vertex(i);
			u.inner_multi_indices(nbrs[i], 0, ind);
			phi[i] = BlockRef(u[ind[0][0]],ind[0][1]);
			if (nonzerophifound==false){
				if (phi[i]!=0){
					nonzerophifound=true;
					nonzerophi=phi[i];
				}
			} else {
				if (nonzerophi*phi[i]<0){
					onls=true;
				}
			}
		};
		if (onls==false) continue;
		// collect neighbor nodes for higher order interpolation
		stageStart[0]=0;
		for (size_t i=1;i<depth+2;i++){
			stageStart[i]=noc;
		}
		//debugUG_LOG("------------------------------------------\n");
		for (size_t stage=0;stage<depth;stage++){
			//debugUG_LOG("stage= " << stage << "\n");
			for (size_t i=stageStart[stage];i<stageStart[stage+1];i++){
				//debugUG_LOG("i =" << i << " stageStart[stage+1]=" << stageStart[stage+1] << "\n");
				CollectNeighbors(nbrCandidates, grid, nbrs[i]);
				for (size_t j=0;j<nbrCandidates.size();j++){
					bool newNeighbor=true;
					for (size_t k=0;k<nbrs.size();k++){
						if (nbrCandidates[j]==nbrs[k]){
							newNeighbor=false;
							break;
						}
					};
					if (newNeighbor==true){
						nbrs.push_back(nbrCandidates[j]);
						//debugUG_LOG("new size : " << nbrs.size() << "\n");
						//debugUG_LOG(nbrs.size()-1 << " " << aaPos[nbrs[nbrs.size()-1]] << "\n");
					};
				};
			};
			stageStart[stage+2]=nbrs.size();
			//debugUG_LOG(" # stageStart[" << stage+2 << "]=" << stageStart[stage+2] << "\n");
		};
		phi.resize(nbrs.size());
		coord.resize(nbrs.size());
		// UG_LOG("nbrs.size=" << nbrs.size() << "\n");
		//debugUG_LOG("### " << aaPos[nbrs[0]] << "\n");
		for (size_t i=0;i<nbrs.size();i++){
			//debug			UG_LOG(i << "\n");
//			UG_LOG("--- " << aaPos[nbrs[i]] << "\n");
			coord[i] = aaPos[nbrs[i]];
			//debugUG_LOG("co0(" << i+1 << ",:)=[" << coord[i][0] << "," << coord[i][1] << "];\n");
			if (i<noc) continue;
			u.inner_multi_indices(nbrs[i], 0, ind);
			phi[i] = DoFRef(u,ind[0]);
			//debug
//			UG_LOG("phi[i]=" << phi[i] << "\n");
		};
		number kappa;
		computeElementCurvature2d(kappa,noc,coord,phi,order);
//		computeElementCurvature2d2(kappa,noc,coord,phi,order,leastSquaresFactor);
		u.multi_indices(elem,1,ind);
		DoFRef(u,ind[0]) = kappa;
		if (m_exactcurvatureknown==true)
			if (abs(kappa+m_exactcurv)>maxnormerr){
				maxnormerr =abs(kappa+m_exactcurv);
			}
		//u.inner_multi_indices(elem, 1, ind);
		//DoFRef(u,ind[1]) = kappa;
	};
	};
	if (m_exactcurvatureknown==true) UG_LOG("curvature maximum error " << maxnormerr << "\n");
	return true;
};

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::computeElementCurvatureFromSides(TGridFunction& u,size_t order,number leastSquaresFactor){
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object elem_type;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator elem_iterator;
	typedef typename elem_type::side side_type;
	typedef typename TGridFunction::template traits<side_type>::const_iterator side_iterator;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& aaPos = u.domain()->position_accessor();
	std::vector<MultiIndex<2> > ind;
	typedef typename Grid::AttachmentAccessor<side_type,ANumber > aSideNumber;
	aSideNumber acEdgeCurvature;
	ANumber aEdgeCurvature;
	grid.template attach_to<side_type>(aEdgeCurvature);
	acEdgeCurvature.access(grid,aEdgeCurvature);
	static number undefined = 1953853528.340591483532;
	SetAttachmentValues(acEdgeCurvature, grid.template begin<side_type>(), grid.template end<side_type>(), undefined);
	number maxnormerr = 0;
	//	loop elements of dimension
	for (int si=0;si<u.num_subsets();++si){
		side_iterator iter = u.template begin<side_type>(si);
		side_iterator iterEnd = u.template end<side_type>(si);
		if (u.num_fct(si)<2){
			UG_THROW("No curvature component in approximation space.");
		}
		if (u.local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, 1)){
			UG_THROW("First component in approximation space must be of Lagrange 1 type.");
		}
		if (u.local_finite_element_id(1) != LFEID(LFEID::PIECEWISE_CONSTANT,0)){
			UG_THROW("Second component in approximation space must be of piecewise constant type.");
		}
		std::vector<MathVector<dim> > coord;
		std::vector<number> phi;
		std::vector<VertexBase*> nbrs;
		std::vector<VertexBase*> nbrCandidates;
		size_t depth=order;
		std::vector<size_t> stageStart(depth+2);
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			side_type* elem = *iter;
			//	get position accessor
			size_t noc=elem->num_vertices();
			nbrs.resize(noc);
			phi.resize(noc);
			// check if zero level set is on element
			bool onls=false;
			bool nonzerophifound=false;
			number nonzerophi;
			for (size_t i=0;i<noc;i++){
				nbrs[i]=elem->vertex(i);
				u.inner_multi_indices(nbrs[i], 0, ind);
				phi[i] = BlockRef(u[ind[0][0]],ind[0][1]);
				if (nonzerophifound==false){
					if (phi[i]!=0){
						nonzerophifound=true;
						nonzerophi=phi[i];
					}
				} else {
					if (nonzerophi*phi[i]<0){
						onls=true;
					}
				}
			};
			if (onls==false) continue;
			// collect neighbor nodes for higher order interpolation
			stageStart[0]=0;
			for (size_t i=1;i<depth+2;i++){
				stageStart[i]=noc;
			}
			for (size_t stage=0;stage<depth;stage++){
				for (size_t i=stageStart[stage];i<stageStart[stage+1];i++){
					CollectNeighbors(nbrCandidates, grid, nbrs[i]);
					for (size_t j=0;j<nbrCandidates.size();j++){
						bool newNeighbor=true;
						for (size_t k=0;k<nbrs.size();k++){
							if (nbrCandidates[j]==nbrs[k]){
								newNeighbor=false;
								break;
							}
						};
						if (newNeighbor==true){
							nbrs.push_back(nbrCandidates[j]);
						};
					};
				};
				stageStart[stage+2]=nbrs.size();
			};
			phi.resize(nbrs.size());
			coord.resize(nbrs.size());
			for (size_t i=0;i<nbrs.size();i++){
				coord[i] = aaPos[nbrs[i]];
				if (i<noc) continue;
				u.inner_multi_indices(nbrs[i], 0, ind);
				phi[i] = DoFRef(u,ind[0]);
			};
			number kappa;
//			if (computeElementCurvature2d2(kappa,noc,coord,phi,order,leastSquaresFactor)==false) UG_THROW("curvature calculation failed.\n");
			computeElementCurvature2d(kappa,noc,coord,phi,order);
			acEdgeCurvature[elem] = kappa;
		}
		elem_iterator elemIter = u.template begin<elem_type>(si);
		elem_iterator elemIterEnd = u.template end<elem_type>(si);
		for(  ;elemIter !=elemIterEnd; ++elemIter)
		{
			//	get Elem
			elem_type* elem = *elemIter;
			typename grid_type::template traits<side_type>::secure_container sides;
			grid.associated_elements(sides, elem );
			bool onls = false;
			number ecurvature = 0;
			size_t nInterSides = 0;
			for (size_t i=0;i<sides.size();i++){
				if (acEdgeCurvature[sides[i]]!= undefined){
					onls = true;
					ecurvature+=acEdgeCurvature[sides[i]];
					nInterSides++;
				}
			}
			if (onls==false) continue;
			number kappa = (number)ecurvature/nInterSides;
			u.multi_indices(elem,1,ind);
			DoFRef(u,ind[0]) = kappa;
			if (m_exactcurvatureknown==true)
				if (abs(kappa+m_exactcurv)>maxnormerr){
					maxnormerr =abs(kappa+m_exactcurv);
				}
		}
	};
	if (m_exactcurvatureknown==true) UG_LOG("curvature maximum error " << maxnormerr << "\n");
	return true;
};
// limit previously computed gradient so that the control-volume-wise linear interpolation function does not introduce new maxima or minima into the data 
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::limit_grad(TGridFunction& uOld, aaGrad& aaGrad)
{
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
	typedef typename TGridFunction::template traits<EdgeBase>::const_iterator EdgeBaseConstIterator;

	// get grid
	typename domain_type::grid_type& grid = *uOld.domain()->grid();

	position_accessor_type& aaPos = uOld.domain()->position_accessor();

	std::vector<MultiIndex<2> > ind;

	//	create Attachment for scv-volume size
	ANumber aMax;
	ANumber aMin;

	//	attach to grid
	grid.attach_to_vertices(aMin);
	grid.attach_to_vertices(aMax);

	//	get attachment accessor to access values
	Grid::VertexAttachmentAccessor<ANumber> aaMin(grid, aMin);
	Grid::VertexAttachmentAccessor<ANumber> aaMax(grid, aMax);
	for (int si=0;si<uOld.num_subsets();++si){
		for(VertexBaseConstIterator iter = uOld.template begin<VertexBase>(si);
				iter !=uOld.template end<VertexBase>(si); ++iter)
				{
		    	    VertexBase* vrt = *iter;
		    	    MathVector<dim> coord;
		    	    coord = aaPos[vrt];
		    	    //	read indices on vertex
		      	    //	get vector holding all indices on the vertex
				    uOld.inner_multi_indices(vrt, 0, ind);
				    aaMax[vrt] = BlockRef(uOld[ind[0][0]],ind[0][1]);
				    aaMin[vrt] = BlockRef(uOld[ind[0][0]],ind[0][1]);
			    }
	}
	for (int si=0;si<uOld.num_subsets();++si){
		//UG_LOG("si " << si << "\n");
		for(EdgeBaseConstIterator iter = uOld.template begin<EdgeBase>(si) ;
				iter !=uOld.template end<EdgeBase>(si); ++iter)
		{
			EdgeBase* edge = *iter;
			VertexBase* vi=edge->vertex(0);
			VertexBase* vj=edge->vertex(1);
			uOld.inner_multi_indices(vi, 0, ind);
			number ui = BlockRef(uOld[ind[0][0]],ind[0][1]);
			uOld.inner_multi_indices(vj, 0, ind);
			number uj = BlockRef(uOld[ind[0][0]],ind[0][1]);
			//UG_LOG("edge " << aaPos[vi] << "-" << aaPos[vj] << " [" << ui << " " << uj << "]\n");
			if (uj<aaMin[vi]){
				aaMin[vi]=uj;
			}
			if (uj>aaMax[vi]){
				aaMax[vi]=uj;
			}
			if (ui<aaMin[vj]){
				aaMin[vj]=ui;
			}
			if (ui>aaMax[vj]){
				aaMax[vj]=ui;
			}
		}

	};
	for (int si=0;si<2;++si){
		for(VertexBaseConstIterator iter = uOld.template begin<VertexBase>(si) ;
				iter !=uOld.template end<VertexBase>(si); ++iter)
		{
/*		    VertexBase* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
			// read indices on vertex
			//	get vector holding all indices on the vertex
			uOld.inner_multi_indices(vrt, 0, ind);
			//UG_LOG(coord << " min=" << aaMin[vrt] << " max=" << aaMax[vrt] << "\n");*/
		}
	}
	for (int si=0;si<1;++si){
	    for(EdgeBaseConstIterator iter = uOld.template begin<EdgeBase>(si) ;
	    		iter !=uOld.template end<EdgeBase>(si); ++iter)
   	    {
         EdgeBase* edge = *iter;
         VertexBase* vi=edge->vertex(0);
         VertexBase* vj=edge->vertex(1);
         MathVector<dim> coordi,coordj,coordij,distVec,gradi,gradj;
         gradi = aaGrad[vi];
         gradj = aaGrad[vj];
		 coordi = aaPos[vi];
		 coordj = aaPos[vj];
		 uOld.inner_multi_indices(vi, 0, ind);
		 number ui = BlockRef(uOld[ind[0][0]],ind[0][1]);
		 uOld.inner_multi_indices(vj, 0, ind);
		 number uj = BlockRef(uOld[ind[0][0]],ind[0][1]);
		 VecScaleAdd(coordij,0.5,coordi,0.5,coordj);
		 VecSubtract(distVec, coordij,coordi);
		 number uij = ui + distVec*gradi;
		 number alpha = 1;
		 if (uij>ui){
		     if (uij>aaMax[vi]) alpha=(aaMax[vi]-ui)/(distVec*gradi);
		     if (alpha<1){
		    	 //UG_LOG("edge " << coordi << " " << coordj << "\n");
		         //UG_LOG(coordi << " u " << ui << " uij "  << uij << " max " << aaMax[vi] << " alpha " << alpha << "\n");
			     aaGrad[vi]*=alpha;
		     };
		 }else{
		     if (uij<aaMin[vi]) alpha=(aaMin[vi]-ui)/(distVec*gradi);
		     if (alpha<1){
		    	 //UG_LOG("edge " << coordi << " " << coordj << "\n");
		    	 //UG_LOG(coordi << " u " << ui << " uij "  << uij << " min " << aaMax[vi] << " alpha " << alpha << "\n");
		    	 aaGrad[vi]*=alpha;
		     };
		 };
		 VecSubtract(distVec, coordij,coordj);
		 uij = uj + distVec*gradj;
		 alpha = 1;
		 if (uij>uj){
		    if (uij>aaMax[vj]) alpha=(aaMax[vj]-uj)/(distVec*gradj);
		    if (alpha<1){
		    	//UG_LOG("-- edge " << coordi << " " << coordj << "\n");
		    	//UG_LOG(coordj << " u " << uj << " uij "  << uij << " max " << aaMax[vj] << " alpha " << alpha << "\n");
		 	    aaGrad[vj]*=alpha;
		    };
		 }else{
		 	if (uij<aaMin[vj]) alpha=(aaMin[vj]-uj)/(distVec*gradj);
		 	if (alpha<1){
		 		//UG_LOG("-- edge " << coordi << " " << coordj << "\n");
		 		//UG_LOG(coordj << " u " << uj << " uij "  << uij << " min " << aaMax[vj] << " alpha " << alpha << "\n");
		 	    aaGrad[vj]*=alpha;
		 	};
		 };
         // UG_LOG(" coord vertex 0 " << aaPos[v0] << " coord vertex 1 " << aaPos[v1] << "\n");
	}

	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	//	coord and vertex array
	MathVector<dim> grad[domain_traits<dim>::MaxNumVerticesOfElem];
	//	get iterators
	ElemIterator iter = uOld.template begin<ElemType>(si);
	ElemIterator iterEnd = uOld.template end<ElemType>(si);
	//	loop elements of dimension
	for(  ;iter !=iterEnd; ++iter)
	{
	    //	get Elem
		ElemType* elem = *iter;

		//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& aaPos = uOld.domain()->position_accessor();

		//	compute center of mass
		MathVector<dim> center;
		center=0;
		size_t noc=elem->num_vertices();
		number u[noc];
    	for(size_t i = 0; i < noc; ++i){
    		vVrt[i] = elem->vertex(i);
    		coCoord[i] = aaPos[vVrt[i]];
			grad[i] = aaGrad[vVrt[i]] ;
			VecAppend(center,coCoord[i]);
			uOld.inner_multi_indices(vVrt[i], 0, ind);
			u[i]=BlockRef(uOld[ind[0][0]],ind[0][1]);
		};
		center/=noc;
		for (size_t i=0;i<noc;++i){
			number alpha=1;
			MathVector<dim> distVec;
			number uCenter;
    		VecSubtract(distVec,center,coCoord[i]);
			uCenter = u[i] + distVec*grad[i];
			if (uCenter>u[i]){
		        if (uCenter>aaMax[vVrt[i]]) alpha=(aaMax[vVrt[i]]-u[i])/(distVec*grad[i]);
		        if (alpha<1){
		        	// UG_LOG("* " << coCoord[i] << " uCenter " << uCenter << " ui " << u[i] << " max " << aaMax[vVrt[i]] << " alpha " << alpha << "\n");
			        aaGrad[vVrt[i]]*=alpha;
		        };
			}else{
		        if (uCenter<aaMin[vVrt[i]]) alpha=(aaMin[vVrt[i]]-u[i])/(distVec*grad[i]);
		        if (alpha<1){
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


// assemble element using upwind
template<typename TGridFunction>
template <typename TElem>
bool FV1LevelSetDisc<TGridFunction>::assemble_element(TElem& elem, DimFV1Geometry<dim>& geo, grid_type& grid,TGridFunction& uNew,const TGridFunction& uOld,aaGrad& aaGradient, aaVol& aaVolume )
{
	//	get domain
	domain_type& domain = *uNew.domain().get();

	//	get grid of domain
	// typename domain_type::grid_type& grid = domain.grid();

	//	create Multiindex
	std::vector<MultiIndex<2> > multInd;

	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

	//	hard code function (fct=0)
	//\todo: generalize
//	size_t fct=0;

	//	id of shape functions used
//	LFEID id = uOld.local_finite_element_id(fct);

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& aaPos = domain.position_accessor();

    //	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get vertices and extract corner coordinates
	const size_t numVertices = elem->num_vertices();
	for(size_t i = 0; i < numVertices; ++i){
		vVrt[i] = elem->vertex(i);
		coCoord[i] = aaPos[vVrt[i]];
	};
	// update fv geometry
	geo.update(elem, &(coCoord[0]), uOld.domain()->subset_handler().get());

    //UG_LOG("geo.num_scvf() " << geo.num_scvf() << "\n");
    //UG_LOG("geo.num_scv() " << geo.num_scv() << "\n");

//  fill node value vector
	std::vector<number> uValue(geo.num_scv());
	size_t noc = geo.num_scv();
	for (size_t i=0;i < noc;i++){
		// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
		uOld.inner_multi_indices(vVrt[i], 0, multInd);
		uValue[i]=BlockRef(uOld[multInd[0][0]],multInd[0][1]);
		//UG_LOG(coCoord[i] << "corner "<< i << " value " << uValue[i] << "\n");
	}

	//UG_LOG("uValues " << uValue[0] << " " << uValue[1] << " " << uValue[2] << "\n");

//  fill grad vector
    MathVector<dim> grad[maxNumCo];
	for (size_t i=0;i < noc;i++){
		// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
		uOld.inner_multi_indices(vVrt[i], 0, multInd);
		grad[i]= aaGradient[vVrt[i]];
		//UG_LOG("corner " << i << " gradient " << grad[i] << "\n");
	}
	
	//  fill corner velocity vector
	MathVector<dim> coVelocity[maxNumCo];
	for (size_t i=0;i<noc;++i){
		coVelocity[i]=0;
	}
	if (m_gamma!=0){
		const int si = 0;
		if (m_imVelocity->requires_grid_fct()){
			//	create storage
			LocalIndices localind;
			LocalVector localu;

			// 	get global indices
			uNew.indices(elem, localind);

			// 	adapt local algebra
			localu.resize(localind);

			// 	read local values of u
			GetLocalVector(localu, uNew);

			//	compute data
	//		try{
				(*m_imVelocity)(coVelocity, geo.scv_global_ips(), m_time, si,
				localu, elem,
				coCoord, geo.scv_local_ips(),
				geo.num_sh());
	//		}
	//		UG_CATCH_THROW("assemble_element : Cannot evaluate velocity data.");
		} else {
			// see user_data.h : 410
			(*m_imVelocity)(coVelocity, geo.scv_global_ips(), m_time, si, geo.num_sh());
		}

        if (m_gamma!=1) for (size_t i=0;i < noc;i++) coVelocity[i]*=m_gamma;
	}
    if (m_delta!=0){
    	for (size_t i=0;i < noc;i++){
    	    number vnorm = VecLength(grad[i]);
    	    if (vnorm>1e-15) for (int j=0;j<dim;j++) coVelocity[i][j] += m_delta/vnorm*grad[i][j];
    	    //UG_LOG("corner " << i << " gradient" << grad[i] << " velocity " << coVelocity[i] << "\n");
    	};
    };
//  fill ipVel velocity vector
	MathVector<dim> ipVelocity[maxNumCo];
	//UG_LOG("num scvf " << geo.num_scvf() << "\n");
//	if ((m_interpolate_v_in_ip==true)||(m_delta!=0)){
		for (size_t ip=0;ip < geo.num_scvf();ip++){
			ipVelocity[ip] = 0;
			const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
			for (size_t co=0;co < noc;co++){
			   // UG_LOG("co " << co << "\n");
				for (int j=0;j<dim;j++){
				    ipVelocity[ip][j] += scvf.shape(co)*coVelocity[co][j];
				    // UG_LOG("ip " << ip << " shape " << co << "[" << j << "]=" << scvf.shape(co) << " co velocity " << j << "=" << coVelocity[co][j]<<"\n");
				};
			};
		}
//	} else
	//{
		//const int si = 0;
		//(*m_imVelocity)(coVelocity, geo.scv_global_ips(), m_time, si, geo.num_sh());;
	//};
    //  fill source vector
	std::vector<number> coSource(noc);
	const int si = 0;
	(*m_imSource)(&coSource[0], geo.scv_global_ips(), m_time, si, geo.num_sh());
    // get finite volume geometry
	// FV1Geometry<TElem,dim> geo;

    //	typename TGridFunction::vector_type& v_vec = *dynamic_cast<typename TGridFunction::vector_type*>(&u);

	size_t base;
	number flux;
	MathVector<dim> distVec;
	// compute fluxes
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	    MathVector<dim> bNode;
	    // 	get current SCVF
	    const typename DimFV1Geometry<dim>::SCVF& scvf = geo.scvf(ip);
	    MathVector<dim>	ipCoord = scvf.global_ip();
	    size_t from = scvf.from();
	    size_t to   = scvf.to();
	    if (scvf.normal()*ipVelocity[ip]>0){
		    base = from;
		} else {
		    base = to;
		};
	    //UG_LOG("from " << coCoord[from] << " to " << coCoord[to] << "\n");
		VecSubtract(distVec, ipCoord,coCoord[base]);
		//UG_LOG("ip vel " << ipVelocity[ip] << " normal " << scvf.normal() << " v*n " << ipVelocity[ip]*scvf.normal() << "\n");
		//UG_LOG("distVec*grad[base] " << distVec*grad[base] << " coSource[base] " << coSource[base] << " grad[base]*coVelocity[base] " << grad[base]*coVelocity[base] << "\n");
		//UG_LOG("u^{n+0.5}_{ip_i} = " << uValue[base] + distVec*grad[base] + 0.5*m_dt*(coSource[base] - grad[base]*coVelocity[base]) << "\n");
		// flux = v * n * u_{ip(i)}^{n+0.5}
		flux = m_dt*(ipVelocity[ip]*scvf.normal())*( uValue[base] + (distVec*grad[base]) + 0.5*m_dt*(coSource[base] - (grad[base]*coVelocity[base])) );
		//UG_LOG(ip << " flux=" << flux << "\n");
		uOld.inner_multi_indices(vVrt[from], 0, multInd);
        BlockRef(uNew[multInd[0][0]],multInd[0][1])-=flux/aaVolume[ vVrt[from] ];
		if (m_divFree==false){
		    BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ];
		    //BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.0*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ];
		};
		//UG_LOG("source flux from " << m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from])))/aaVolume[ vVrt[from] ] << "\n");
		//UG_LOG("v*n=" << (ipVelocity[ip]*scvf.normal()) << " uExtra=" <<  (uValue[from] + 0.5*m_dt*(coSource[from] - (grad[from]*coVelocity[from]))) << " vol=" << aaVolume[ vVrt[from] ] << "\n");
		uOld.inner_multi_indices(vVrt[to], 0, multInd);
        BlockRef(uNew[multInd[0][0]],multInd[0][1])+=flux/aaVolume[ vVrt[to] ];
		if (m_divFree==false){
		    BlockRef(uNew[multInd[0][0]],multInd[0][1])-=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to])))/aaVolume[ vVrt[to] ];
			//BlockRef(uNew[multInd[0][0]],multInd[0][1])-=m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[to] + 0.0*m_dt*(coSource[to] - (grad[to]*coVelocity[to])))/aaVolume[ vVrt[to] ];
		};
		//UG_LOG("source flux to " << m_dt*(ipVelocity[ip]*scvf.normal())*(uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to])))/aaVolume[ vVrt[to] ] << "\n");
		//UG_LOG("v*n=" << (ipVelocity[ip]*scvf.normal()) << " uExtra=" <<  (uValue[to] + 0.5*m_dt*(coSource[to] - (grad[to]*coVelocity[to]))) << " vol=" << aaVolume[ vVrt[to] ] << "\n");
        number localCFL = std::max(m_dt*abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[from] ],m_dt*abs(ipVelocity[ip]*scvf.normal())/aaVolume[ vVrt[to] ] );
        //UG_LOG("localCFL " << localCFL << "\n");
        if (localCFL>m_maxCFL){
            m_maxCFL = localCFL;
        };
        //UG_LOG("global CFL " << m_maxCFL << "\n");
         //		number dist = VecTwoNorm(distVec);
	};
	// boundary
	//	evaluate finite volume geometry

		//UG_LOG("nr of bdry subsets: " << geo.num_boundary_subsets() << " nr of bdry faces: " << geo.num_bf() << "\n");
		if (geo.num_bf()>0){
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
				    for (size_t co=0;co<noc;co++){
					    //UG_LOG("corner " << co << " num_sh " << bf.num_sh() << "\n");
				        for (int j=0;j<dim;j++){
				            bipVelocity[j] += bf.shape(co)*coVelocity[co][j];
		                    bipGrad[j] += bf.shape(co) * globalGradVec[co][j];
				    	};
				        bipU += bf.shape(co) *uValue[co];
				        bipSource += bf.shape(co) * coSource[co];
				    }
    				//flux = m_dt*(bipVelocity*bf.normal())*( bipU + 0.5*m_dt*(bipSource - (bipGrad*bipVelocity)) );
    				flux = m_dt*(bipVelocity*bf.normal())*uValue[nodeID];// first order approximation
    				uOld.inner_multi_indices(vVrt[nodeID], 0, multInd);
    				DoFRef(uNew, multInd[0])-=flux/aaVolume[ vVrt[nodeID] ];
    				//UG_LOG("coord=" << bf.global_ip( )<< "vel=" << bipVelocity << " n=" << bf.normal() << " flux=" << flux << " source flux=" << m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] + 0.5*m_dt*(coSource[nodeID] - (grad[nodeID]*coVelocity[nodeID])))/aaVolume[ vVrt[nodeID] ] << "\n");
    				//UG_LOG("v*n=" << (bipVelocity*bf.normal()) << " uExtra=" << (uValue[nodeID] + 0.5*m_dt*(coSource[nodeID] - (grad[nodeID]*coVelocity[nodeID]))) << " volume=" << aaVolume[ vVrt[nodeID] ] << "\n");
    				if (m_divFree==false){
    				    //BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] + 0.5*m_dt*(coSource[nodeID] - (grad[nodeID]*coVelocity[nodeID])))/aaVolume[ vVrt[nodeID] ];
    				    DoFRef(uNew, multInd[0])+=m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] )/aaVolume[ vVrt[nodeID] ];// first order approximation
    				    //BlockRef(uNew[multInd[0][0]],multInd[0][1])+=m_dt*(bipVelocity*bf.normal())*(uValue[nodeID] + 0.5*m_dt*(bipSource - (grad[nodeID]*coVelocity[nodeID])))/aaVolume[ vVrt[nodeID] ];
    				};
		    	};
		     };
		 };

	// give out corner values for debug
//	for (size_t i=0;i < noc;i++){
			// if (dd.template inner_multi_indices<VertexBase>(vVrt[i], 0, multInd) != 1) return false;
//			dd.inner_multi_indices(vVrt[i], 0, multInd);
//			uValue[i]=BlockRef(uNew[multInd[0][0]],multInd[0][1]);
			//UG_LOG(coCoord[i] << "NEW corner "<< i << " value " << uValue[i] << "\n");
//	}

	//UG_LOG("NEW uValues " << uValue[0] << " " << uValue[1] << " " << uValue[2] << "\n");
	// end of "for debug"

//	we're done
	return true;
}

/*************************

COMPUTE VOLUME OF CONTROL VOLUMES

**************************/
template <typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::
calculate_vertex_vol(TGridFunction& u,aaVol& aaVolume)
{
	//	get domain
		domain_type& domain = *u.domain().get();

	//	create a FV Geometry for the dimension
		DimFV1Geometry<dim> geo;

	//	element iterator type
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
		typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

	//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
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
				};
			}
		}

	return true ;

}

// compute gradient in vertices and volume of control volume
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::
calculate_vertex_grad_vol(TGridFunction& u, aaGrad& aaGradient,aaVol& aaVolume)
{
	//	domain type
	//	typedef typename TGridFunction::domain_type domain_type;

	/// dof distribution type
	//	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;

	// 	dimension of world (and reference element)
	//	static const int dim = domain_type::dim;

	//	get domain
		domain_type& domain = *u.domain().get();

	//	get grid of domain
		typename domain_type::grid_type& grid = *domain.grid();

	//	create Multiindex
		std::vector<MultiIndex<2> > multInd;

	//	create a FV Geometry for the dimension
		DimFV1Geometry<dim> geo;

	//	get element iterator type
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
		typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

	//	hard code function (fct=0)
	//\todo: generalize
		size_t fct=0;

		// initialize attachment value
		SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
		SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

		//	coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	sum up all contributions of the sub control volumes to one vertex in an attachment
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
	        if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
	        if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)==true) continue;
	        if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)==true) continue;
		//	get iterators
			ElemIterator iter = u.template begin<ElemType>(si);
			ElemIterator iterEnd = u.template end<ElemType>(si);

		//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
			//	get Elem
				ElemType* elem = *iter;

			//	get position accessor
				typedef typename domain_type::position_accessor_type position_accessor_type;
				const position_accessor_type& aaPos = domain.position_accessor();

			//	get vertices and extract corner coordinates
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = aaPos[vVrt[i]];
				};

			//	evaluate finite volume geometry
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

				//UG_LOG("Num Verts loaded: "<<vVrt.size()<<"\n");
				//UG_LOG("Num SCV computed: "<<geo.num_scv()<<"\n");m_vPP.

				number uValue[domain_traits<dim>::MaxNumVerticesOfElem];

			//	read indices on vertex
				size_t noc = geo.num_scv();
				for (size_t i=0;i < noc;i++)
				{
				//	get indices of function fct on vertex
					u.inner_multi_indices(vVrt[i], fct, multInd);

				//	read value of index from vector
					uValue[i]=BlockRef(u[multInd[0][0]],multInd[0][1]);

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
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
	position_accessor_type aaPos = u.domain()->position_accessor();
	for (int si=0;si < u.num_subsets();++si){
		//UG_LOG("si " << si << "\n");
	    for(VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
						   iter != u.template end<VertexBase>(si); ++iter)
	    {
	    //	get vertex
		    VertexBase* vrt = *iter;
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
	return true ;
}

// compute gradient in vertices and volume of control volume in region given by sign of level set function
template<typename TGridFunction>
	bool FV1LevelSetDisc<TGridFunction>::
	calculate_vertex_grad_vol_sign(TGridFunction& u, aaGrad& aaGradient,aaVol& aaVolume,TGridFunction& phi,int sign)
	{
		//	domain type
		//	typedef typename TGridFunction::domain_type domain_type;

		/// dof distribution type
		//	typedef typename TGridFunction::dof_distribution_type dof_distribution_type;

		// 	dimension of world (and reference element)
		//	static const int dim = domain_type::dim;

		//	get domain
		domain_type& domain = *u.domain().get();

		//	get grid of domain
		typename domain_type::grid_type& grid = *domain.grid();

		//	create Multiindex
		std::vector<MultiIndex<2> > multInd;

		//	create a FV Geometry for the dimension
		DimFV1Geometry<dim> geo;

		//	get element iterator type
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
		typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

		//	hard code function (fct=0)
		//\todo: generalize
		size_t fct=0;

		// initialize attachment value
		SetAttachmentValues(aaVolume, grid.vertices_begin(), grid.vertices_end(), 0);
		SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0);

		//	coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//	sum up all contributions of the sub control volumes to one vertex in an attachment
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
	        if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
	        if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)==true) continue;
	        if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)==true) continue;
			//	get iterators
			ElemIterator iter = u.template begin<ElemType>(si);
			ElemIterator iterEnd = u.template end<ElemType>(si);

			//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
				//	get Elem
				ElemType* elem = *iter;

				//	get position accessor
				typedef typename domain_type::position_accessor_type position_accessor_type;
				const position_accessor_type& aaPos = domain.position_accessor();

				//	get vertices and extract corner coordinates
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = aaPos[vVrt[i]];
				};

				//	evaluate finite volume geometry
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

				//UG_LOG("Num Verts loaded: "<<vVrt.size()<<"\n");
				//UG_LOG("Num SCV computed: "<<geo.num_scv()<<"\n");m_vPP.

				number uValue[domain_traits<dim>::MaxNumVerticesOfElem];

				size_t noc = geo.num_scv();
				bool rightsign=true;
				for (size_t i=0;i < noc;i++)
				{
					//	get indices of function fct on vertex
					u.inner_multi_indices(vVrt[i], fct, multInd);

					//	read value of index from vector
					uValue[i]=BlockRef(u[multInd[0][0]],multInd[0][1]);
					if (sign==-1){
						if (BlockRef(phi[multInd[0][0]],multInd[0][1])>0){
							rightsign = false;
							break;
						}
				    };
					if (sign==1){
						if (BlockRef(phi[multInd[0][0]],multInd[0][1])>0){
							rightsign = false;
							break;
						};
					};
					//	debug log
					//UG_LOG("corner " << i << " " << uValue[i] << "\n");
				}
				if (rightsign==false) continue;
				for (size_t i=0;i < noc;i++)
				{
					//	get indices of function fct on vertex
					u.inner_multi_indices(vVrt[i], fct, multInd);

					//	read value of index from vector
					uValue[i]=BlockRef(u[multInd[0][0]],multInd[0][1]);

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
		typedef typename domain_type::position_accessor_type position_accessor_type;
		typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
		position_accessor_type aaPos = u.domain()->position_accessor();
		for (int si=0;si < u.num_subsets();++si){
			//UG_LOG("si " << si << "\n");
			for(VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
				iter != u.template end<VertexBase>(si); ++iter)
			{
				//	get vertex
				VertexBase* vrt = *iter;
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
		return true ;
	}

// set Dirichlet values in solution vector
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::assign_dirichlet(TGridFunction& numsol){
	//	get domain of grid function
		domain_type& domain = *numsol.domain().get();

	//	get grid type of domain
		typedef typename domain_type::grid_type grid_type;

		typedef typename domain_type::position_accessor_type position_accessor_type;

	//	UG_LOG("dirichlet\n");

		// UG_LOG("nr dir ss " << m_dirichlet_sg.size() << "\n");

		typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
		for(size_t i = 0; i < m_dirichlet_sg.size(); ++i)
		{
	        const int si = m_dirichlet_sg[i];
	        // UG_LOG("Dirichlet boundary is: "<<si<< "\n");
			for(VertexBaseConstIterator iter = numsol.template begin<VertexBase>(si);
										   iter != numsol.template end<VertexBase>(si); ++iter)
			{
			//	get vertex
				VertexBase* vrt = *iter;
				number exactVal;
				position_accessor_type aaPos = domain.position_accessor();

			//	get vector holding all indices on the vertex
				std::vector<MultiIndex<2> > ind;

				const size_t numInd = numsol.inner_multi_indices(vrt, 0, ind);

			//	check indices
				if(numInd != 1) {UG_LOG("ERROR: Wrong number of indices!"); return false;}

				(*m_imDirichlet)(&exactVal,&aaPos[vrt],m_time,si,1);
				BlockRef(numsol[ind[0][0]],ind[0][1]) = exactVal;

				// MathVector<dim> coord = aaPos[vrt];
				// UG_LOG("coord " << coord[0] << "," << coord[1] << " <> " << BlockRef(numsol[ind[0][0]],ind[0][1]) << "\n");

				//if ((coord[0]==-1)||(coord[0]==1)||(coord[1]==-1)||(coord[1]==1)){
			    	//BlockRef(numsol[ind[0][0]],ind[0][1]) = exactVal;
				//};
				//if ((coord[0]==0)||(coord[0]==1)||(coord[1]==0)||(coord[1]==1)){
				//};
		     }
		};
	return true;
}

// in region given by level set sign:
// overwrite values in parameter unew with values of uold
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::overwrite(TGridFunction& unew,TGridFunction& uold,TGridFunction& phi,int sign)
{
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

	for (int si=0;si<unew.num_subsets();++si){
		for(VertexBaseConstIterator iter = unew.template begin<VertexBase>(si);
									   iter != unew.template end<VertexBase>(si); ++iter)
		{
			VertexBase* vrt = * iter;

		//	read indices on vertex
			std::vector<MultiIndex<2> > ind;
			unew.inner_multi_indices(vrt, 0, ind);
		    number phiValue = BlockRef(phi[ind[0][0]],ind[0][1]);
			int nodeSign;
			if (phiValue<0)  nodeSign =-1;
			if (phiValue>0)  nodeSign = 1;
			if (phiValue==0) nodeSign = 0;
		    if (nodeSign==sign)	BlockRef(unew[ind[0][0]],ind[0][1]) = BlockRef(uold[ind[0][0]],ind[0][1]);
		}
	};
    return true;
};

// in region given by level set sign:
// overwrite values in parameter unew with parameter value
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::overwrite(TGridFunction& unew,number value,TGridFunction& phi,int sign)
{
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;

	for (int si=0;si<unew.num_subsets();++si){

		for(VertexBaseConstIterator iter = unew.template begin<VertexBase>(si);
									   iter != unew.template end<VertexBase>(si); ++iter)
		{
			VertexBase* vrt = * iter;

		//	read indices on vertex
			std::vector<MultiIndex<2> > ind;
			unew.inner_multi_indices(vrt, 0, ind);
		    number phiValue = BlockRef(phi[ind[0][0]],ind[0][1]);
			int nodeSign;
			if (phiValue<0)  nodeSign =-1;
			if (phiValue>0)  nodeSign = 1;
			if (phiValue==0) nodeSign = 0;
		    if (nodeSign==sign)	BlockRef(unew[ind[0][0]],ind[0][1]) = value;
		}
	};
    return true;
};

// for test problems:
// compute error by using analytical solution
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::compute_error(TGridFunction& numsol)
{
//	get domain of grid function
	domain_type& domain = *numsol.domain().get();

//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

//	get grid of domain
	grid_type& grid = *domain.grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

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
	if (calculate_vertex_vol(numsol,aaVolume)==false){UG_LOG("ERROR: gradient computation failed in compute_error function!"); };

	//UG_LOG("----------------------------\n");

	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
	if(!bRes) {UG_LOG("Error while calculating CV Volume.\n"); return false;}
	for (int si=0;si<numsol.num_subsets();++si){
		// UG_LOG("*** " << si << "\n");
		for(VertexBaseConstIterator iter = numsol.template begin<VertexBase>(si);
									   iter != numsol.template end<VertexBase>(si); ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			number exactVal;
			position_accessor_type aaPos = domain.position_accessor();

		//	get vector holding all indices on the vertex
			std::vector<MultiIndex<2> > ind;

			const size_t numInd = numsol.multi_indices(vrt, 0, ind);

		//	check indices
			if(numInd != 1) {UG_LOG("ERROR: Wrong number of indices!"); return false;}

			(*m_imDirichlet)(&exactVal,&aaPos[vrt],m_time,si,1);
			number differ = abs(BlockRef(numsol[ind[0][0]],ind[0][1])-exactVal);
		
			l1Error += aaVolume[vrt] * differ;
			l2Error += aaVolume[vrt] * differ*differ;

			if (m_print==true){
			 //   if (differ>0)
			    UG_LOG("coord=" << aaPos[vrt] << " value=" << BlockRef(numsol[ind[0][0]],ind[0][1]) << " exact=" << exactVal << " error=" << differ << "\n");
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

// initialize level set function with analytical solution
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::init_function(TGridFunction& u)
{
//	get domain of grid function
	domain_type& domain = *u.domain().get();

//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	typedef typename domain_type::position_accessor_type position_accessor_type;

	//	read indices on vertex
	position_accessor_type aaPos = domain.position_accessor();

	typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
	for (int si=0;si<domain.subset_handler()->num_subsets();++si){
		for(VertexBaseConstIterator iter = u.template begin<VertexBase>(si);
									   iter != u.template end<VertexBase>(si); ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			MathVector<dim> coord;
			coord = aaPos[vrt];
		//	get vector holding all indices on the vertex
			std::vector<MultiIndex<2> > ind;
			u.inner_multi_indices(vrt, 0, ind);
			BlockRef(u[ind[0][0]],ind[0][1]) = analytic_solution(m_time,coord);
	     }
	};
	return true;
};

// for runtime testing, delete later
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::runtimetest(TGridFunction& uNew)
{
	//	get domain of grid function
	domain_type& domain = *uNew.domain().get();

	//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	//	get grid of domain
//	grid_type& grid = domain.grid();

	//	element iterator type
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;
	
 	//	create a FV Geometry for the dimension
	DimFV1Geometry<dim> geo;

//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& aaPos = domain.position_accessor();

//	resize corners
//	std::vector<MathVector<dim> > coCoord;
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	for (int si=0;si<uNew.num_subsets();++si){
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
	for (int si=0;si<uNew.num_subsets();++si){
		//	get iterators
		ElemIterator iter = uNew.template begin<ElemType>(si);
		ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;
			//	get vertices of the Elem
			std::vector<VertexBase*> vVrt;
			CollectVertices(vVrt, grid, elem);
			//	get position accessor
			typedef typename domain_type::position_accessor_type position_accessor_type;
			const position_accessor_type& aaPos = domain.position_accessor();
			//	resize corners
			std::vector<MathVector<dim> > coCoord;
			//	extract corner coordinates
			for(size_t i = 0; i < vVrt.size(); ++i)
				coCoord.push_back( aaPos[vVrt[i]] );
		};
	};
	

	for (int si=0;si<domain.subset_handler().num_subsets();++si){
		//	get iterators
		ElemIterator iter = uNew.template begin<ElemType>(si);
		ElemIterator iterEnd = uNew.template end<ElemType>(si);
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			ElemType* elem = *iter;
			//	get vertices of the Elem
			std::vector<VertexBase*> vVrt;
			CollectVertices(vVrt, grid, elem);
		};
	};

	MathVector<dim> coord;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	position_accessor_type aaPos = domain.position_accessor();

	for (int si=0;si<uNew.num_subsets();++si){
		VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
	    VertexBaseConstIterator iterEnd = uNew.template end<VertexBase>(si);
		for(;iter != iterEnd; ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;
			coord = aaPos[vrt];
		}
	};*/
	return true;
}

// assign subsets as given by level set function 
// ug subset system would have to be changed for this to be useful
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::update_ls_subsets(TGridFunction& phi){
	//	get domain of grid function
    domain_type& domain = *phi.domain().get();

    //	get grid type of domain
    typedef typename domain_type::grid_type grid_type;

    //	get grid of domain
    // grid_type& grid = domain.grid();

    typedef typename domain_type::position_accessor_type position_accessor_type;

	//	element iterator type
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

    //	create Multiindex
    std::vector<MultiIndex<2> > ind;

	//	get element iterator type
	m_inactive_sg.set_subset_handler(domain.subset_handler());
    for (int si=0;si<domain.subset_handler()->num_subsets();++si){
    	UG_LOG("******************* si " << si << " **********************\n");
		if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
        if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)==true) continue;
		ElemIterator iter = phi.template begin<ElemType>(si);
     	ElemIterator iterEnd = phi.template end<ElemType>(si);
     	//	loop elements of dimension
     	for(  ;iter !=iterEnd; ++iter)
     	{
  	    //	get Elem
        	ElemType* elem = *iter;

        	//	coord and vertex array
        	VertexBase* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

        	//	get vertices and extract corner coordinates
        	const size_t numVertices = elem->num_vertices();
        	for(size_t i = 0; i < numVertices; ++i){
        		vVrt[i] = elem->vertex(i);
        	};
            //	get position accessor
      		typedef typename domain_type::position_accessor_type position_accessor_type;

	        //	resize corners
	        std::vector<MathVector<dim> > coCoord;
     		//	compute center of mass
     		MathVector<dim> center;
     		std::vector<MathVector<dim> > grad;
     		center=0;
      		int noc=elem->num_vertices();
     		number phiCo[noc];
           	for(int i = 0; i < noc; ++i){
	        	phi.inner_multi_indices(vVrt[i], 0, ind);
     			phiCo[i]=BlockRef(phi[ind[0][0]],ind[0][1]);
      		};
			int firstNonzero=-1;
			for (int j=0;j<noc;j++){
			    if (phiCo[j]!=0){
				    firstNonzero=j;
				    break;
				};
			};
			if (firstNonzero==-1) // all element nodes are on zero ls
			     domain.subset_handler()->assign_subset(elem,m_inside_elements_si);
			// add nodes on zero ls to zero ls node subsetgroup
			for (int i=0;i<noc;i++){
				if (phiCo[i]==0){
    			    int oldindex = domain.subset_handler()->get_subset_index(vVrt[i]);
	            	if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(oldindex)==true) continue;
		            if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(oldindex)==true) continue;
                    domain.subset_handler()->assign_subset(vVrt[i],m_onls_nodes_si);
 				};
			};
			bool onls=false;
			for (int i=firstNonzero+1;i<noc;i++){
			    if (phiCo[firstNonzero]*phiCo[i]<0){
			    	onls=true;
			    	break;
			    }
			};
			if (onls==false){
			    if (phiCo[firstNonzero]<0){
			    //	UG_LOG("si before " << domain.subset_handler().get_subset_index(elem));
				    domain.subset_handler()->assign_subset(elem,m_inside_elements_si);
				    UG_LOG("element is inside \n");
				//    UG_LOG("si after " << domain.subset_handler().get_subset_index(elem) << "\n");
				//    UG_LOG("-- nr of subsets: " << phi.num_subsets() << " " << m_inside_elements_si <<  "\n");
				};
				if (phiCo[firstNonzero]>0){
				//	UG_LOG("si before " << domain.subset_handler().get_subset_index(elem));
				    domain.subset_handler()->assign_subset(elem,m_outside_elements_si);
				//    UG_LOG("si after " << domain.subset_handler().get_subset_index(elem) <<  "\n");
				//    UG_LOG("|| nr of subsets: " << phi.num_subsets() << " " << m_outside_elements_si << "\n");
				};
			};
			if (onls==true){
				//UG_LOG("si before " << domain.subset_handler().get_subset_index(elem));
			    domain.subset_handler()->assign_subset(elem,m_onls_elements_si);
			    //UG_LOG("si after " << domain.subset_handler().get_subset_index(elem) <<  "\n");
			    //UG_LOG("// nr of subsets: " << phi.num_subsets() << " " << m_onls_elements_si << "\n");
			    for (int i=0;i<noc;i++){
   					int oldindex = domain.subset_handler()->get_subset_index(vVrt[i]);
	    		    if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(oldindex)==true) continue;
       	            if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(oldindex)==true) continue;
				    if (phiCo[i]<0){
					    domain.subset_handler()->assign_subset(vVrt[i],m_inside_nodes_si);
					    UG_LOG("node is inside \n");
					};
					if (phiCo[i]>0){
					    domain.subset_handler()->assign_subset(vVrt[i],m_outside_nodes_si);
					};
					if (phiCo[i]==0){
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
       	 for(VertexBaseConstIterator iter = phi.template begin<VertexBase>(sindex);
    	   iter != phi.template end<VertexBase>(sindex); ++iter)
       	 {
       	     VertexBase* vrt = *iter;
       	    ++count;
       	 }
    	 UG_LOG(count << " nodes in subset\n");
    };*/
  	return true;
}

template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::init_ls_subsets(TGridFunction& phi){
    create_ls_subsets(phi);
	if (update_ls_subsets(phi)==false) return false;
    return true;
}

template<typename TGridFunction>
void FV1LevelSetDisc<TGridFunction>::create_ls_subsets(TGridFunction& phi){
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

// main function for assembling and solving ls equation as described in 
// Frolkovic/Mikula, HIGH-RESOLUTION FLUX-BASED LEVEL SET METHOD, SIAM
template<typename TGridFunction>
bool FV1LevelSetDisc<TGridFunction>::advect_lsf(TGridFunction& uNew,TGridFunction& uOld)
{
	//	get domain of grid function
	domain_type& domain = *uNew.domain().get();

	//	get grid type of domain
	typedef typename domain_type::grid_type grid_type;

	//	get grid of domain
	grid_type& grid = *domain.grid();

	typedef typename domain_type::position_accessor_type position_accessor_type;

	//	element iterator type
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::geometric_base_object ElemType;

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


	// const dof_distribution_type& dd = uNew.dof_distribution();
	m_maxCFL=0;
	
	//UG_LOG("***************************************************\n");
    //UG_LOG("***************************************************\n");
	//UG_LOG("***************************************************\n");

	/*for (int si=0;si<uNew.num_subsets();++si){
		for(VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
								   iter != grid.template end<VertexBase>(si); ++iter)
	    {
		    VertexBase* vrt = *iter;
		    typedef typename TGridFunction::std::vector<MultiIndex<2> > index_type;

	   //	get vector holding all indices on the vertex
		    index_type ind;

    	//	read indices on vertex
		    if (dd.inner_multi_indices(vrt, 0, ind)!=1){UG_LOG("ERROR: Wrong number of indices!"); return false;};

		    BlockRef(uNew[ind[0][0]],ind[0][1]) = BlockRef(uOld[ind[0][0]],ind[0][1]);
		    //UG_LOG("uNew: " << BlockRef(uNew[ind[0][0]],ind[0][1]) << "uOld: " << BlockRef(uOld[ind[0][0]],ind[0][1]) << "\n");
	    }
	};*/
	VecAssign(uNew,uOld);
    MathVector<dim> coord;
    position_accessor_type aaPos = domain.position_accessor();

    //	get vector holding all indices on the vertex
	std::vector<MultiIndex<2> > ind;
    //	read indices on vertex
	for (size_t step=0;step<m_nrOfSteps;step++)
	{
	    // calculate scv size and gradient
	    if (calculate_vertex_grad_vol(uNew,aaGradient, aaVolume)==false){UG_LOG("ERROR: gradient computation failed!"); };
	    // SetAttachmentValues(aaGradient, grid.vertices_begin(), grid.vertices_end(), 0); for debug set gradient to 0
	    if (m_limiter==true){
	    	limit_grad(uNew,aaGradient);
	    };
	    // UG_LOG("num_subsets: " << uOld.num_subsets() << "\n");
	    for (int si=0;si<uOld.num_subsets();++si){
	        //UG_LOG("si " << si << "\n");
	        if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
	        if (m_neumann_sg.size()!=0) if (m_neumann_sg.contains(si)==true) continue;
	        if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)==true) continue;
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

		    //	loop elements of dimension
		    for(  ;iter !=iterEnd; ++iter)
		    {
		        //	get Elem
			    //UG_LOG("*** ELEM ***\n");
			    ElemType* elem = *iter;
			    //UG_LOG("element \n");
			    // uNew = uOld
			    assemble_element(elem, geo, *uNew.domain()->grid(),uNew,uOld,aaGradient,aaVolume);
		    };
	    };
		typedef typename TGridFunction::template traits<VertexBase>::const_iterator VertexBaseConstIterator;
		for (int si=0;si<uOld.num_subsets();++si){
			if (m_dirichlet_sg.size()!=0) if (m_dirichlet_sg.contains(si)==true) continue;
			if (m_inactive_sg.size()!=0) if (m_inactive_sg.contains(si)==true) continue;
			std::vector<number> sourceValue(1);
			MathVector<dim>* sourceCo;
			for(VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
					    			                      iter != uNew.template end<VertexBase>(si); ++iter){
				VertexBase* vrt = *iter;
				sourceCo[0]= aaPos[vrt];
				(*m_imSource)(&sourceValue[0],sourceCo,m_time,si,1);
				uNew.inner_multi_indices(vrt, 0, ind);
				BlockRef(uNew[ind[0][0]],ind[0][1]) += m_dt*sourceValue[0];
			};
	    }
	    m_time += m_dt;
	    UG_LOG("m_dt " << m_dt << " m_time " << m_time << "\n");
	    m_timestep_nr++;
	    assign_dirichlet(uNew);
	    // overwrite inactive nodes with old solution
	    for(size_t i = 0; i < m_inactive_sg.size(); ++i)
	    {
	    	 const int si = m_inactive_sg[i];
	    	 UG_LOG("inactive si: " << si << "\n");
	    	 for(VertexBaseConstIterator iter = uNew.template begin<VertexBase>(si);
	    	 									   iter != uNew.template end<VertexBase>(si); ++iter)
	    	 {
	    	     VertexBase* vrt = *iter;
	    	     UG_LOG("*\n");
	    	     uNew.inner_multi_indices(vrt, 0, ind);
	    		 BlockRef(uNew[ind[0][0]],ind[0][1]) = BlockRef(uOld[ind[0][0]],ind[0][1]);
	    	 }
	    };

	    UG_LOG("step: " << m_timestep_nr << "\n");
	    UG_LOG("time: " << m_time << "\n");
        UG_LOG("max CFL nr " << m_maxCFL << "\n");
	    //};
	    if (m_nrOfSteps>1){
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

} // end namespace LevelSet
} // end namespace ug

#endif /* LEVEL_SET_UTIL_IMPL_H_ */
