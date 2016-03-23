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

Matrix::Matrix(){

}

Matrix::Matrix(unsigned v, unsigned h, unsigned min, unsigned max){
    this->vsize = v;
    this->hsize  = h;
    this->M = new double*[this->vsize];
    for(int i=0; i<this->vsize;i++){
        this->M[i]=new double[this->hsize];
    }
    for(unsigned i=0;i<v;i++)
        for(unsigned j=0;j<h;j++)
            M[i][j]=r_random(min, max);
}

double r_random(double min, double max){
    min = min+0.1;
    max = max-0.1;
    std::random_device rd;
    std::uniform_int_distribution<double> uid(min, max);
    return uid(rd);
}

int Matrix::Get_vsize() const{
    return this->vsize;
}

int Matrix::Get_hsize() const{
    return this->hsize;
}

void Matrix::Add(int i, int j, double value){
    if(i<=vsize && j<=hsize){
        this->M[i][j]=value;
    }else{
        std::cout<<"Size error";
    }
}

double Matrix::Get(int i, int j) const{
    return this->M[i][j];
}

void Matrix::Show() const{
    for(int i=0;i<this->vsize;i++){
        for(int j=0;j<this->hsize;j++){
            std::cout << this->M[i][j] << "\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void Matrix::SwapRows(int k, int i){
    double tmp;
    for(int j=0;j<this->hsize;j++){
        tmp=M[i][j];
        M[i][j]=M[k][j];
        M[k][j]=tmp;
    }
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
        std::cout<<"Multipliing"<<std::endl;
        
        for(unsigned i = 0; i < aVsize; i++)
            for(unsigned j = 0; j < bHsize; j++)
                for(unsigned k = 0; k < aHsize; k++){
                    double tmp = res.Get(i, j)+Get(i, k) * B.Get(k,j);
                    res.Add(i, j, tmp);
                }
        
        
    }
    return res;
}
#endif /* Matrix_cpp */