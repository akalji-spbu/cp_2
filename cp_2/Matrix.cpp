//
//  Matrix.cpp
//  cp_2
//
//  Created by Nikolay Tikhonov on 22.03.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#ifndef Matrix_cpp
#define Matrix_cpp

#include "Matrix.h"

Matrix::Matrix(std::ifstream& infile){
    infile >> this->vsize;
    infile >> this->hsize;
    this->M = new double*[this->vsize];
    for(int i=0; i<this->vsize;i++){
        this->M[i]=new double[this->hsize];
    }
    for(int i=0;i<this->vsize;i++){
        for(int j=0;j<this->hsize;j++){
            infile >> this->M[i][j];
        }
    }
}

Matrix::Matrix(Matrix const &A){
    vsize = A.Get_vsize();
    hsize = A.Get_hsize();
    this->M = new double*[this->vsize];
    for(int i=0; i<this->vsize;i++){
        this->M[i]=new double[this->hsize];
    }
    for(int i=0;i<this->vsize;i++){
        for(int j=0;j<this->hsize;j++){
            this->M[i][j]=A.Get(i,j);
        }
    }
}


Matrix::Matrix(int vsize, int hsize){
    this->vsize = vsize;
    this->hsize  = hsize;
    this->M = new double*[this->vsize];
    for(int i=0; i<this->vsize;i++){
        this->M[i]=new double[this->hsize];
    }
    for(int i=0;i<this->vsize;i++){
        for(int j=0;j<this->hsize;j++){
            this->M[i][j]=0.0;
        }
    }
}

Matrix::Matrix(unsigned N){
    this->vsize = N;
    this->hsize  = N;
    this->M = new double*[N];
    for(int i=0; i<N;i++){
        this->M[i]=new double[N];
    }
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            this->M[i][j]=0.0;
        }
    }
}

Matrix::~Matrix(){
    if (M != nullptr) {
        for (int i = 0; i < vsize; ++i)
        {
            if (M[i] != nullptr) delete[] M[i];
        }
        delete[] M;
    }
}

Matrix::Matrix(){

}

int Matrix::Get_vsize() const{
    return this->vsize;
}

int Matrix::Get_hsize() const{
    return this->hsize;
}

void Matrix::Add(unsigned i, unsigned j, double value){
    if(i<=vsize && j<=hsize){
        this->M[i][j]=value;
    }else{
        std::cout<<"Size error";
    }
}


Matrix Matrix::operator- (const Matrix& A) const{
    if (this->Get_vsize() != A.Get_vsize() || this->Get_hsize() != A.Get_hsize()) throw;
    Matrix B(this->Get_vsize(),this->Get_hsize());
    
    for (int i = 0; i < this->Get_vsize(); ++i)
        for (int j = 0; j < this->Get_hsize(); ++j)
            B.Add(i,j, this->Get(i,j) - A.Get(i,j));
    
    return B;
}

Matrix Matrix::operator+ (const Matrix& A) const{
    if (this->Get_vsize() != A.Get_vsize() || this->Get_hsize() != A.Get_hsize()) throw;
    Matrix B(this->Get_vsize(),this->Get_hsize());
    
    for (int i = 0; i < this->Get_vsize(); ++i)
        for (int j = 0; j < this->Get_hsize(); ++j)
            B.Add(i,j, this->Get(i,j)+A.Get(i,j));
    
    return B;
}

Matrix Matrix::operator*(const Matrix &B) const {
    
    unsigned aVsize = Get_vsize();
    unsigned aHsize = Get_hsize();
    unsigned bVsize = B.Get_vsize();
    unsigned bHsize = B.Get_hsize();
    
    unsigned rVsize, rHsize;
    rVsize=aVsize;
    rHsize=bHsize;
    
    Matrix res(rVsize, rHsize);
    
    for(unsigned i=0;i<rVsize;i++){
        for(unsigned j=0;j<rHsize;j++){
            res.Add(i, j, 0.0);
        }
    }
    
    if(aHsize==bVsize){
        for(unsigned i = 0; i < aVsize; i++)
            for(unsigned j = 0; j < bHsize; j++)
                for(unsigned k = 0; k < aHsize; k++){
                    double tmp = res.Get(i, j)+Get(i, k) * B.Get(k,j);
                    res.Add(i, j, tmp);
                }
    }
    return res;
}


double Matrix::Get(int i, int j) const{
    return this->M[i][j];
}



void Matrix::Show() const{
    for(int i=0;i<this->vsize;i++){
        for(int j=0;j<this->hsize;j++){
            std::cout << this->M[i][j] << "\t \t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void Matrix::SwapColumns(unsigned a, unsigned b){
    unsigned N = this->Get_vsize();
    if (a>N || b>N) throw;
    double tmp;
    for (int i = 0; i < N; ++i){
        tmp = M[i][a];
        M[i][a] = M[i][b];
        M[i][b] = tmp;
    }
    return;
}

void Matrix::SwapRows(unsigned a, unsigned b){
    if (a == b) return;
    double * tmp = M[a];
    M[a] = M[b];
    M[b] = tmp;
    return;
}

Matrix& Matrix::operator= (const Matrix& A){
    
    M = new double*[A.vsize];
    vsize = A.vsize;
    hsize = A.hsize;
    
    for (int i = 0; i < A.vsize; ++i)
        M[i] = new double[A.hsize];
    
    for (int i = 0; i < A.vsize; ++i)
        for (int j = 0; j < A.hsize; ++j)
            M[i][j] = A.M[i][j];
    
    return *this;
}


void DiagonalizeMatrix(Matrix &M){
    unsigned N = M.Get_vsize();
        for(unsigned i=0;i<N;i++){
            double tmp=0.0;
            for(unsigned j=0;j<N;j++)
                tmp+=fabs(M.Get(i, j));
            if(M.Get(i, i)<0) M.Add(i, i, -tmp);
            else M.Add(i, i, tmp);
        }
}

void RandomizeMatrix(Matrix &M, double min, double max){
    unsigned n = M.Get_vsize();
    unsigned m = M.Get_hsize();
    for(unsigned i=0;i<n;i++)
        for(unsigned j=0;j<m;j++)
            M.Add(i, j, r_random(min,max));
}

double r_random(double min, double max){
    std::random_device rd;
    std::uniform_real_distribution<double> uid(min, max);
    return uid(rd);
}

Matrix transpose(Matrix &M){
    unsigned n = M.Get_vsize();
    unsigned m = M.Get_hsize();
    Matrix A(M.Get_hsize(),M.Get_vsize());
    for(unsigned i=0;i<n;i++)
        for (unsigned j=0; j<m; j++)
            A.Add(j,i,M.Get(i,j));
    return A;
}

void Matrix::clear(){
    for(unsigned i=0;i<vsize;i++)
        for(unsigned j=0;j<hsize;j++)
            M[i][j]=0.0;
}

void Matrix::insertDiag(double A){
    for(unsigned i = 0; i<vsize; ++i)
        for(unsigned j = 0; j<hsize; ++j){
            if(i==j) M[i][j]=A;
    }
}

double Matrix::norm(){
    unsigned n=vsize;
    unsigned m=hsize;
    double max = 0;
    for (unsigned i=0; i<n; ++i){
        double tmp=0;
        for (int j=0; j<m; ++j)
            tmp+=fabs(M[i][j]);
        if (tmp>max)max=tmp;
    }
    return max;
}



#endif /* Matrix_cpp */