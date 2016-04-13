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
#include <random>

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
    
    ~Matrix();
    
    
    Matrix operator * (const Matrix& B) const;
    Matrix operator - (const Matrix& A) const;
    Matrix operator + (const Matrix& A) const;
    
    Matrix& operator= (const Matrix& A);
    void Add(unsigned i, unsigned j, double value);
    int Get_vsize() const;
    int Get_hsize() const;
    void Show() const;
    double Get(int i, int j) const;
    Matrix transpose() const;
    void SwapRows(unsigned a, unsigned b);
    void SwapColumns(unsigned a, unsigned b);
    void clear();
    void insertDiag(double A);
    double norm();
};

double r_random(double min, double max);
void DiagonalizeMatrix(Matrix &M);
void RandomizeMatrix(Matrix &M, double min, double max);
Matrix transpose(Matrix &M);


#endif /* Matrix_hpp */
