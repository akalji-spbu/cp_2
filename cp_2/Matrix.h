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
    Matrix(unsigned v, unsigned h, unsigned min, unsigned max);
    
    
    Matrix operator * (const Matrix& B) const;
    
    void Add(int i, int j, double value);
    int Get_vsize() const;
    int Get_hsize() const;
    void Show() const;
    double Get(int i, int j) const;
    void SwapRows(int p, int i);
};

double r_random(double min, double max);


#endif /* Matrix_hpp */
