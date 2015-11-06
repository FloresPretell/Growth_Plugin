/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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

namespace ug{
namespace LevelSet{

int bubblesort(std::vector<number> array,std::vector<int>& list)
{
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
bool multMatVec(const std::vector<number>& avec, const std::vector<number>& b,std::vector<number>& c,size_t m,size_t n)
{
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
            const std::vector<number>& b /* rhs */)
{
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
		max=std::abs(A[i][i]);
		for (k=i+1;k<n;k++){
			if (std::abs(A[k][i])>max){
				j=k;
				max=std::abs(A[k][i]);
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
	i = n;
	do {
		i--;
		number s=z[i];
		for (j=n-1;j>i;j--){
			s=s-x[j]*U[i][j];
		};
		x[i]=s/U[i][i];
	} while (i != 0);
	return true;
};


template<size_t n>
bool solveLS(MathVector<n>& x,const MathMatrix<n,n>& M,const MathVector<n>& b)
{
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
bool leastSquares(MathVector<n>& x,const MathMatrix<m,n>& M,const MathVector<m>& b)
{
	if(m<n){
		UG_THROW("Least squares method not suitable for m x n matrix with m<n.\n");
	}
	MathMatrix<n,n> M2;
	MathVector<n> b2;
	MatMultiplyMTM(M2,M);
	TransposedMatVecMult(b2,M,b);
	return solveLS(x,M2,b2);
}

bool leastSquares(std::vector<number>& x,const std::vector<number>& mField,const std::vector<number>& b)
{
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


// averages positions by arithmetic mean
/*
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
bool LevelSetCurvature<TGridFunction>::computeElementCurvature2d(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order)
{
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
    }
	else if (order==3){
		dxphi=coeffs[1]+2*coeffs[3]*elemBasePoint[0]+coeffs[4]*elemBasePoint[1]+3*coeffs[6]*elemBasePoint[0]*elemBasePoint[0]+2*coeffs[7]*elemBasePoint[0]*elemBasePoint[1]+coeffs[8]*elemBasePoint[1]*elemBasePoint[1];
		dyphi=coeffs[2]+coeffs[4]*elemBasePoint[0]+2*coeffs[5]*elemBasePoint[1]+coeffs[7]*elemBasePoint[0]*elemBasePoint[0]+2*coeffs[8]*elemBasePoint[0]*elemBasePoint[1]+3*coeffs[9]*elemBasePoint[1]*elemBasePoint[1];
    } else {
    	UG_THROW("Case unsupported");
    }
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
        if (std::abs(phival)<1e-12) break;
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
bool LevelSetCurvature<TGridFunction>::computeElementCurvature2d2(number& kappa,size_t elementnoc,const std::vector<MathVector<dim> >& co,std::vector<number> phi,size_t order,number nodefactor)
{
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
        if (std::abs(phival)<1e-12) break;
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
bool LevelSetCurvature<TGridFunction>::computeElementCurvatureOnGrid2d(TGridFunction& u,size_t order,number leastSquaresFactor)
{
	typedef typename TGridFunction::domain_type domain_type;
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object ElemType;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& aaPos = u.domain()->position_accessor();
	std::vector<DoFIndex> ind;
	number maxnormerr = 0;
	//	loop elements of dimension
	for (int si=0;si<u.num_subsets();++si){
		ElemIterator iter = u.template begin<ElemType>(si);
		ElemIterator iterEnd = u.template end<ElemType>(si);
		if (u.num_fct(si)<2){
			UG_THROW("No curvature component in approximation space.");
		}
		if (u.local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, dim, 1)){
			UG_THROW("First component in approximation space must be of Lagrange 1 type.");
		}
		if (u.local_finite_element_id(1) != LFEID(LFEID::PIECEWISE_CONSTANT,dim,0)){
			UG_THROW("Second component in approximation space must be of piecewise constant type.");
		}
		std::vector<MathVector<dim> > coord;
		std::vector<number> phi;
		std::vector<Vertex*> nbrs;
		std::vector<Vertex*> nbrCandidates;
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
			u.inner_dof_indices(nbrs[i], 0, ind);
			phi[i] = DoFRef(u, ind[0]);
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
			u.inner_dof_indices(nbrs[i], 0, ind);
			phi[i] = DoFRef(u, ind[0]);
			//debug
//			UG_LOG("phi[i]=" << phi[i] << "\n");
		};
		number kappa;
		computeElementCurvature2d(kappa,noc,coord,phi,order);
//		computeElementCurvature2d2(kappa,noc,coord,phi,order,leastSquaresFactor);
		u.dof_indices(elem,1,ind);
		DoFRef(u,ind[0]) = kappa;
		if (m_exactcurvatureknown==true)
			if (std::abs(kappa+m_exactcurv)>maxnormerr){
				maxnormerr =std::abs(kappa+m_exactcurv);
			}
		//u.inner_dof_indices(elem, 1, ind);
		//DoFRef(u,ind[1]) = kappa;
	};
	};
	if (m_exactcurvatureknown==true) UG_LOG("curvature maximum error " << maxnormerr << "\n");
	return true;
};

template<typename TGridFunction>
bool LevelSetCurvature<TGridFunction>::computeElementCurvatureFromSides(TGridFunction& u,size_t order,number leastSquaresFactor)
{
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	// get grid
	typename domain_type::grid_type& grid = *u.domain()->grid();
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
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
		if (u.local_finite_element_id(0) != LFEID(LFEID::LAGRANGE, dim, 1)){
			UG_THROW("First component in approximation space must be of Lagrange 1 type.");
		}
		if (u.local_finite_element_id(1) != LFEID(LFEID::PIECEWISE_CONSTANT,dim,0)){
			UG_THROW("Second component in approximation space must be of piecewise constant type.");
		}
		std::vector<MathVector<dim> > coord;
		std::vector<number> phi;
		std::vector<Vertex*> nbrs;
		std::vector<Vertex*> nbrCandidates;
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
				u.inner_dof_indices(nbrs[i], 0, ind);
				phi[i] = DoFRef(u, ind[0]);
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
				u.inner_dof_indices(nbrs[i], 0, ind);
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
			u.dof_indices(elem,1,ind);
			DoFRef(u,ind[0]) = kappa;
			if (m_exactcurvatureknown==true)
				if (std::abs(kappa+m_exactcurv)>maxnormerr){
					maxnormerr =std::abs(kappa+m_exactcurv);
				}
		}
	};
	if (m_exactcurvatureknown==true) UG_LOG("curvature maximum error " << maxnormerr << "\n");
	return true;
};

} // end namespace LevelSet
} // end namespace ug

/* End of File */
