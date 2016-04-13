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

void lu();
void qr();
void ja();
void se();

char filenameA[] = "/Users/akalji/Downloads/A.txt";
char filenameB[] = "/Users/akalji/Downloads/b.txt";

int main(int argc, const char * argv[]) {
    double min,max;
    int hsize,vsize;
    std::cout<<"Please enter vsize of matrix: "; std::cin>> vsize;
    std::cout<<"Please enter hsize of matrix: "; std::cin>> hsize;
    std::cout<<"Please enter min: "; std::cin>> min;
    std::cout<<"Please enter max: "; std::cin>> max;
    Matrix M(vsize,hsize);
    Matrix b(vsize,1);
    RandomizeMatrix(b, min, max);
    RandomizeMatrix(M, min, max);
    DiagonalizeMatrix(M);
    Matrix x(Seidel(M, b));
    std::cout<<"x: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    x.Show();
    std::cout<<"Test: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    Matrix test = M*x-b;
    test.Show();
    Matrix xj(Jacobi(M, b));
    std::cout<<"x: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    xj.Show();
    std::cout<<"Test: "<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    Matrix testj = M*xj-b;
    testj.Show();

    
    return 0;
}

void lu(){
    std::ifstream infile(filenameA);
    Matrix A(infile);
    infile.close();
    Matrix A2(A);
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
    std::ifstream infileB(filenameB);
    Matrix B(infileB);
    std::cout<<"Matrix B"<<std::endl;
    std::cout<<"______________________________________________"<<std::endl;
    B.Show();
    Matrix x(A2.Get_vsize(),1);
    Matrix nB(A2.Get_vsize(),1);
    infile.close();
    
    nB = P1*B;
    SOLE(L,U,nB,x,rank);
    
    
    
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
}

void qr(){
    std::cout<<"1) Read from file"<<std::endl;
    std::cout<<"2) Random diagonal matrix"<<std::endl;
    int mode=2;
    //std::cout<<"Select mode: "; std::cin>> mode;
    switch(mode){
        case 1:{
            std::ifstream infile(filenameA);
            Matrix M(infile);
            infile.close();
            std::ifstream infileB(filenameB);
            Matrix b(infileB);
            infileB.close();
            Matrix M2(M);
            infile.close();
            std::cout<<"Source matrix is: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            M2.Show();
            std::cout<<"b matrix is: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            b.Show();
            Matrix R(M);
            Matrix Q(M.Get_vsize(),M.Get_hsize());
            QR(M,Q,R);
            std::cout<<"Q: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            Q.Show();
            std::cout<<"R: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            R.Show();
            Matrix x(QRSLE(Q, R, b));
            Matrix QtrQ = Q*transpose(Q);
            std::cout<<"Q*transpose Q: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            QtrQ.Show();
            Matrix rtQR = transpose(Q)*R;
            std::cout<<"transpose Q * R: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            rtQR.Show();
            std::cout<<"x: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            x.Show();
            std::cout<<"Test: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            Matrix test = M*x-b;
            test.Show();
            break;
        }
        case 2:{
            double min,max;
            int hsize,vsize;
            std::cout<<"Please enter vsize of matrix: "; std::cin>> vsize;
            std::cout<<"Please enter hsize of matrix: "; std::cin>> hsize;
            std::cout<<"Please enter min: "; std::cin>> min;
            std::cout<<"Please enter max: "; std::cin>> max;
            Matrix M(vsize,hsize);
            Matrix b(vsize,1);
            RandomizeMatrix(b, min, max);
            RandomizeMatrix(M, min, max);
            DiagonalizeMatrix(M);
            Matrix M2(M);
            std::cout<<"Source matrix is: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            M2.Show();
            std::cout<<"b matrix is: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            b.Show();
            Matrix R(M);
            Matrix Q(M.Get_vsize(),M.Get_hsize());
            QR(M,Q,R);
            std::cout<<"Q: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            Q.Show();
            std::cout<<"R: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            R.Show();
            Matrix x(QRSLE(Q, R, b));
            Matrix QtrQ = Q*transpose(Q);
            std::cout<<"Q*transpose Q: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            QtrQ.Show();
            Matrix rtQR = transpose(Q)*R;
            std::cout<<"transpose Q * R: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            rtQR.Show();
            std::cout<<"x: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            x.Show();
            std::cout<<"Test: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            Matrix test = M*x-b;
            test.Show();
            break;
        }
        default:
            break;
    }
}

void ja(){
            double min,max;
            int hsize,vsize;
            std::cout<<"Please enter vsize of matrix: "; std::cin>> vsize;
            std::cout<<"Please enter hsize of matrix: "; std::cin>> hsize;
            std::cout<<"Please enter min: "; std::cin>> min;
            std::cout<<"Please enter max: "; std::cin>> max;
            Matrix M(vsize,hsize);
            Matrix b(vsize,1);
            RandomizeMatrix(b, min, max);
            RandomizeMatrix(M, min, max);
            DiagonalizeMatrix(M);
            Matrix M2(M);
            std::cout<<"Source matrix is: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            M2.Show();
            std::cout<<"b matrix is: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            b.Show();
            Matrix x(Jacobi(M, b));
            std::cout<<"x: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            x.Show();
            std::cout<<"Test: "<<std::endl;
            std::cout<<"______________________________________________"<<std::endl;
            Matrix test = M*x-b;
            test.Show();
}

void se(){
        double min,max;
        int hsize,vsize;
        std::cout<<"Please enter vsize of matrix: "; std::cin>> vsize;
        std::cout<<"Please enter hsize of matrix: "; std::cin>> hsize;
        std::cout<<"Please enter min: "; std::cin>> min;
        std::cout<<"Please enter max: "; std::cin>> max;
        Matrix M(vsize,hsize);
        Matrix b(vsize,1);
        RandomizeMatrix(b, min, max);
        RandomizeMatrix(M, min, max);
        DiagonalizeMatrix(M);
        Matrix M2(M);
        std::cout<<"Source matrix is: "<<std::endl;
        std::cout<<"______________________________________________"<<std::endl;
        M2.Show();
        std::cout<<"b matrix is: "<<std::endl;
        std::cout<<"______________________________________________"<<std::endl;
        b.Show();
        Matrix x(Seidel(M, b));
        std::cout<<"x: "<<std::endl;
        std::cout<<"______________________________________________"<<std::endl;
        x.Show();
        std::cout<<"Test: "<<std::endl;
        std::cout<<"______________________________________________"<<std::endl;
        Matrix test = M*x-b;
        test.Show();
}