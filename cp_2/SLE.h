//
//  SLE.hpp
//  cp_2
//
//  Created by Nikolay Tikhonov on 22.03.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#ifndef SLE_hpp
#define SLE_hpp

#include <iostream>
#include <cmath>
#include "Matrix.h"

void P1P2LU(Matrix &A, Matrix &P1, Matrix &P2, Matrix &L, Matrix &U, unsigned &rank, unsigned &swaps);
unsigned PLU(Matrix &A, Matrix &P,Matrix &L, Matrix &U);
double det(Matrix &U, unsigned SwapsNum);
bool SOLE(Matrix &L, Matrix &U, Matrix &b, Matrix &x, int rank);
void inverse(Matrix& L, Matrix& U, int rank, Matrix &inverse);
double cond(Matrix& A, Matrix& inverse);

#endif /* SLE_hpp */
