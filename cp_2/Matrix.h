//
//  Matrix.hpp
//  cp_2
//
//  Created by Nikolay Tikhonov on 22.03.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//


#ifndef Matrix_hpp
#define Matrix_hpp


#include <iostream>
#include <fstream>
#include <cmath>

class Matrix{
private:
    int vsize, hsize;
    double **M;
public:
    Matrix(std::ifstream& infile);
    Matrix(int vsize, int hsize);
    Matrix(Matrix const &A);
    Matrix();
    Matrix(unsigned N);
    
    Matrix operator * (const Matrix& B) const;
    
    void Add(int i, int j, double value);
    int Get_vsize();
    int Get_hsize();
    void Show();
    double Get(int i, int j);
    void SwapRows(int p, int i);
};

Matrix multiply(Matrix &A, Matrix &B);


#endif /* Matrix_hpp */
