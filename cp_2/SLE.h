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

unsigned PLU(Matrix &A, Matrix &P,Matrix &L, Matrix &U);
double det(Matrix &U, unsigned SwapsNum);

#endif /* SLE_hpp */
