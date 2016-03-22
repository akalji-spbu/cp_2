//
//  CLE.cpp
//  cp_2
//
//  Created by Nikolay Tikhonov on 21.03.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//


#ifndef SLE_cpp
#define SLE_cpp
#include "SLE.h"


unsigned PLU(Matrix &A, Matrix &P,Matrix &L, Matrix &U){
    unsigned n=A.Get_vsize();
    unsigned m=A.Get_hsize();
    unsigned SwapsNum;
    // make identity matrix P
    for(unsigned i=0;i<n;i++){
        for(unsigned j=i;j<m;j++){
            if(j==i)
                P.Add(i,j,1);
            else
                P.Add(i,j,0);
        }
    }
    for( int i = 0; i < n; i++ ){
        double pivotValue = 0;
        int pivot = -1;
        for(int row = i; row < n; row++){
            if( fabs(A.Get(row,i)) > pivotValue ){
                pivotValue = fabs(A.Get(row,i));
                pivot = row;
            }
        }
        P.SwapRows(pivot, i);
        A.SwapRows(pivot, i);
        SwapsNum++;
        double tmp1,tmp2;
        for(unsigned j = i+1; j < m; j++){
            tmp1 = A.Get(j,i)/A.Get(i,i);
            A.Add(j,i,tmp1);
            for(unsigned k = i+1; k < m; k++){
                tmp2 = A.Get(j,k);
                tmp2 -= A.Get(j,i) * A.Get(i,k);
                A.Add(j,k,tmp2);
            }
        }
    }
    for(unsigned i=0;i<n;i++){
        for(unsigned j=i;j<m;j++){
            if(j==i)
                L.Add(i,j,1);
            else
                L.Add(j,i,A.Get(j,i));

            U.Add(i,j,A.Get(i,j));
        }
    }
    return SwapsNum;
}

double det(Matrix &U, unsigned SwapsNum){
    double det = 1;
    if(U.Get_vsize()!=U.Get_hsize()){
        std::cout<<"Is not square matrix"<<std::endl;
        exit(0);
    }
    unsigned N = U.Get_vsize();
    for(int i=0;i<N;i++){
        det*=U.Get(i,i);
    }
    
    if (SwapsNum % 2 == 0)
        return det;
    else
        return -det;
    
}


#endif /* SLE_cpp */