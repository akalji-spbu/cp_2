//
//  main.cpp
//  cp_2
//
//  Created by Nikolay Tikhonov on 22.03.16.
//  Copyright © 2016 Nikolay Tikhonov. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "SLE.h"

int main(int argc, const char * argv[]) {
    std::ifstream infile("/Users/akalji/Desktop/matrix1.txt");
    Matrix A(infile);
    A.Show();
    Matrix P(A.Get_vsize(),A.Get_hsize());
    Matrix L(A.Get_vsize(),A.Get_hsize());
    Matrix U(A.Get_vsize(),A.Get_hsize());
    
    Matrix R;
    Matrix T;
    
    unsigned Swaps=PLU(A, P, L, U);
    
    double det1 = det(U, Swaps);


    P.Show();
    L.Show();
    U.Show();
    
    R = multiply(P, L);
    T = multiply(R, U);
    T.Show();
    std::cout<<"Det A = "<<det1<<std::endl;
    
    return 0;
}
