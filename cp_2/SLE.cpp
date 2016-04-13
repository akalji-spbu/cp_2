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


//PLU METHODS
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void P1P2LU(Matrix &A, Matrix &P, Matrix &Q, Matrix &L, Matrix &U, unsigned &rank, unsigned &swaps){
    swaps = 0;
    unsigned n = A.Get_vsize();
    unsigned m = A.Get_hsize();
    rank = n;
    
    for (unsigned i = 0; i < n; ++i)
        for (unsigned j = 0; j < m; ++j)
            if(i==j){
                P.Add(i,j,1);
                Q.Add(i,j,1);
            }else{
                P.Add(i,j,0);
                Q.Add(i,j,0);
            }
    
    for (unsigned i = 0; i < n; ++i){
        double max = fabs(A.Get(i,i));
        unsigned max_row = i, max_col = i;
        for (unsigned j = i; j < n; ++j)
            for (unsigned k = i; k < n; ++k)
                if (fabs(A.Get(j,k)) > max){
                    max = fabs(A.Get(j,k));
                    max_row = j;
                    max_col = k;
                }
        if (fabs(max)<pow(10,-15)) {
            --rank;
            continue;
        }
        if (i != max_row) {
            A.SwapRows(i, max_row);
            P.SwapRows(i, max_row);
            swaps++;
        }
        if (i != max_col){
            A.SwapColumns(i, max_col);
            Q.SwapColumns(i, max_col);
            swaps++;
        }
        
        for (unsigned j = i + 1; j < n; ++j)
            A.Add(j,i,A.Get(j,i)/A.Get(i,i));
        
        for (unsigned j = i + 1; j < n; ++j)
            for (int k = i + 1; k < n; ++k)
                A.Add(j,k,A.Get(j,k)-A.Get(j,i) * A.Get(i,k));
    }
    
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
        {
            if (i == j)
                L.Add(j,i,1);
            else
                L.Add(j,i,A.Get(j,i));
            U.Add(i,j,A.Get(i,j));
        }
}



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
        if(pivot!=i){
            P.SwapRows(pivot, i);
            A.SwapRows(pivot, i);
            SwapsNum++;
        }
        double tmp1,tmp2;
        for(unsigned j = i+1; j < m; j++){
            if(A.Get(i,i)!=0.0)
                tmp1 = A.Get(j,i)/A.Get(i,i);
            else continue;
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

void SOLE(Matrix &L,Matrix &U,Matrix &b,Matrix &x,int rank){
    unsigned n=L.Get_vsize();
    Matrix y(n,1);
    for (unsigned i=0; i<n;++i){
        double tmp=b.Get(i,0);
        for (unsigned j=0; j<i;++j)
            tmp = tmp-L.Get(i,j)*y.Get(j,0);
        y.Add(i,0,tmp/L.Get(i,i));
    }
    if (rank==n){
        for (int i = n-1; i >= 0; --i){
            double tmp = y.Get(i,0);
            for (unsigned j = n-1; j>i; --j)
                tmp=tmp-U.Get(i,j)*x.Get(j,0);
            x.Add(i,0,tmp/U.Get(i,i));
        }
    }
    else{
        bool res=false;
        for (unsigned i=n-1;i>=rank;--i)
            if (fabs(y.Get(i,0))>=pow(10,-15)) res = true;
        if(res){
            std::cout<<"SOLE incompatible"<<std::endl;
            exit(0);
        }
        else {
            for (int i=n-1;i>=rank;--i) x.Add(i,0,1);
            for (int i=rank-1; i>=0;--i){
                double tmp=y.Get(i,0);
                for (int j=n-1; j>i;--j)
                    tmp=tmp-U.Get(i,j)*x.Get(j,0);
                x.Add(i,0,tmp/U.Get(i,i));
            }
        }
    }
}

void inverse(Matrix &L, Matrix &U, int rank, Matrix &inverse){
    int n = L.Get_vsize();
    Matrix *x, b(n,1), tmp(n,1);
    x = new Matrix[n];
    for (int i=0; i<n; ++i){
        b.Add(i,0,1);
        if (i>0) b.Add(i-1,0, 0);
        SOLE(L,U,b,tmp,rank);
        x[i]=tmp;
        x[i].Show();
    }
    
    //for (int i = 0; i < n; ++i){
        //std::cout<<"Xi test______________________\n";
       // x[i].Show();}

    
    //std::cout<<"Get test\n";
    for (int i = 0; i < n; ++i){
        //std::cout<<"Xi test\n";
        x[i].Show();
        for (int j = 0; j < n; ++j){
            
            //std::cout<<x[i].Get(j,0)<<std::endl;
        
            inverse.Add(j,i,x[i].Get(j,0));}
    }}

double cond(Matrix& A, Matrix& inverse){
    int  n = A.Get_vsize();
    double num1 = 0, num2 = 0;
    double tmp = 0;
    for (int i = 0; i < n; ++i){
        tmp = 0;
        for (int j = 0; j < n; ++j)
            tmp += fabs(A.Get(i,j));
        if (num1 < tmp) num1 = tmp;
    }
    tmp = 0;
    for (int i = 0; i < n; ++i){
        tmp = 0;
        for (int j = 0; j < n; ++j)
            tmp += fabs(inverse.Get(i,j));
        if (num2 < tmp) num2 = tmp;
    }
    return num1*num2;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//END PLU METHODS


//QR METHODS
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void QR(Matrix &A, Matrix &Q, Matrix &R){
    unsigned n = A.Get_vsize();
    unsigned m = A.Get_hsize();
    Q.insertDiag(1);
    for(int i = 0; i < n; i++){
        for(int j = i+1; j<m; j++){
            Matrix tmpQ(n, m);
            tmpQ.insertDiag(1);
            if (fabs(R.Get(j,i)) < 1e-15) continue;
            double cos = R.Get(i,i) / sqrt(pow(R.Get(i,i), 2) + pow(R.Get(j,i), 2));
            double sin = R.Get(j,i) / sqrt(pow(R.Get(i,i), 2) + pow(R.Get(j,i), 2));
            tmpQ.Add(i,i,cos);
            tmpQ.Add(i,j,sin);
            tmpQ.Add(j,i,-sin);
            tmpQ.Add(j,j,cos);
            Q = tmpQ*Q;
            R = tmpQ*R;
        }
    }
}

Matrix QRSLE(Matrix &Q, Matrix &R, Matrix &b){
    unsigned n = Q.Get_vsize();
    unsigned m = Q.Get_vsize();
    Matrix x(n,1);
    Matrix y(n,1);
    y = Q*b;
    int rank = n-1;
    while (fabs(R.Get(rank,rank)) < 1e-13) --rank;
    for (int i = n - 1; i > rank; --i)
        x.Add(i,0,1);
    for (int i = rank; i >= 0; --i){
        x.Add(i,0,y.Get(i,0));
        for (int j = i + 1; j < m; ++j)
            x.Add(i,0,x.Get(i,0) - x.Get(j,0) * R.Get(i,j));
        x.Add(i,0,x.Get(i,0) /R.Get(i,i));
    }
    return x;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//END QR METHODS

//SEIDEL METHODS
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix Seidel(Matrix &M, Matrix &b){
    unsigned n = M.Get_vsize();
    unsigned m = M.Get_hsize();
    Matrix TMP(n, 1), B(n), c(n, 1), x(n, 1);
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < m; ++j)
            if (i == j) {
                c.Add(i,0,b.Get(i,0)/M.Get(i,i));
                B.Add(i,i,0);
            }
            else B.Add(i,j,-M.Get(i,j)/M.Get(i,i));
    }
    double r=0;
    for (int i=0; i<n; ++i){
        double tmp = 0;
        for (int j=i+1; j<m; ++j)
            tmp += fabs(B.Get(i,j));
        if(tmp>r) r=tmp;
    }
    
    double q=B.norm();
    TMP=c;
    unsigned k=0;
    while (true){
        
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < n; ++j)
                if(j < i) x.Add(i,0, x.Get(i,0)+ B.Get(i,j) * x.Get(j,0));
                else x.Add(i,0,(x.Get(i,0)+B.Get(i,j) * TMP.Get(j,0)));
            x.Add(i,0,(x.Get(i,0) + c.Get(i,0)));
        }
        k++;
        if (fabs(r*(x-TMP).norm()/(1-q))<1e-13) break;
        TMP = x;
        for (int i=0; i<n; ++i)
            x.Add(i,0,0);
    }
    std::cout<<"Iterations of Seidel: "<<k<<std::endl;
    return x;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//END SEIDEL METHODS


//JACOBI METHODS
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix Jacobi(Matrix &M, Matrix&b){
    unsigned n = M.Get_vsize();
    unsigned m = M.Get_hsize();
    Matrix TMP(n, 1), B(n,n), c(n, 1), x(n, 1); //Матрицы заплнены нулями
    for (int i=0; i<n; ++i)
        for (int j=0; j<m; ++j)
            if(i==j){
                c.Add(i,0,b.Get(i,0)/M.Get(i,i));
                B.Add(i,i,0);
            }
            else B.Add(i,j,(-M.Get(i,j)/M.Get(i,i)));
    TMP = c;
    double q = B.norm();
    unsigned k=0;
    while(true){
        x = B*TMP+c;
        double tmp = q*(x-TMP).norm()/(1-q);
        k++;
        if (fabs(tmp) <= 1e-13) break;
        TMP = x;
    }
    std::cout<<"Iterations of Jacobi: "<<k<<std::endl;
    return x;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//END JACOBI METHODS

#endif /* SLE_cpp */