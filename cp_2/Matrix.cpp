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

int Matrix::Get_vsize(){
    return this->vsize;
}

int Matrix::Get_hsize(){
    return this->hsize;
}

void Matrix::Add(int i, int j, double value){
    if(i<=vsize && j<=hsize){
        this->M[i][j]=value;
    }else{
        std::cout<<"Size error";
    }
}

double Matrix::Get(int i, int j){
    return this->M[i][j];
}

void Matrix::Show(){
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

Matrix multiply(Matrix &A, Matrix &B){
    
    unsigned aVsize = A.Get_vsize();
    unsigned aHsize = A.Get_hsize();
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
                    double tmp = res.Get(i, j)+A.Get(i, k) * B.Get(k,j);
                    res.Add(i, j, tmp);
                }
        
        
    }
    return res;
}



#endif /* Matrix_cpp */