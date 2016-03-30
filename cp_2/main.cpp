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
    std::ifstream infile("/Users/akalji/Downloads/A.txt");
    Matrix A(infile);
    infile.close();
    std::ifstream infileA("/Users/akalji/Downloads/A.txt");
    Matrix A2(infileA);
    infile.close();
    std::cout<<"Source matrix is: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    A2.Show();
    Matrix P1(A.Get_vsize(),A.Get_hsize());
    Matrix P2(A.Get_vsize(),A.Get_hsize());
    Matrix L(A.Get_vsize(),A.Get_hsize());
    Matrix U(A.Get_vsize(),A.Get_hsize());
    unsigned rank, swaps;
    std::cout<<"Let's start PLU decopposition: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    
    P1P2LU(A, P1, P2, L, U, rank, swaps);
    std::cout<<"P1"<<std::endl;
    P1.Show();
    std::cout<<"P2"<<std::endl;
    P2.Show();
    std::cout<<"L"<<std::endl;
    L.Show();
    std::cout<<"U"<<std::endl;
    U.Show();
    
    Matrix P1AP2;
    Matrix LU;
    
    
    P1AP2=P1*A2*P2;
    std::cout<<"P1*A*P2"<<std::endl;
    P1AP2.Show();
    LU=L*U;
    std::cout<<"L*U"<<std::endl;
    LU.Show();
    
    std::cout<<"Det A is: "<< det(U, swaps) <<std::endl;
    std::cout<<"Rank A is: "<< rank <<std::endl;
    
    std::cout<<"Let's start SOLE: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    std::ifstream infileB("/Users/akalji/Downloads/b.txt");
    Matrix B(infileB);
    std::cout<<"Matrix B"<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    B.Show();
    Matrix x(A2.Get_vsize(),1);
    Matrix nB(A2.Get_vsize(),1);
    infile.close();
    
    nB = P1*B;
    if (!SOLE(L,U,nB,x,rank)) {
        std::cout << "SOLE incompatible" << std::endl;
        return 0;
    }
    x=P2*x;
    std::cout<<"X: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    x.Show();
    
    std::cout<<"test: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    Matrix Re;
    Re=A2*x-B;
    Re.Show();
    
    std::cout<<"Let's start inverse: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    Matrix hinverse(A2.Get_vsize());
    inverse(L,U,rank,hinverse);
    Matrix inverse = P2*hinverse*P1;
    inverse.Show();
    std::cout<<"inverse*A: "<<std::endl;
    Matrix uninverse = A2*inverse;
    uninverse.Show();
    
    std::cout<<"Let's start calculate conditional number: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    double condit;
    condit = cond(A2, inverse);
    std::cout << condit<<std::endl;
    
    
    return 0;
}
