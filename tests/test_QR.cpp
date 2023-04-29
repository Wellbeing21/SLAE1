#include <gtest/gtest.h>
#include "../Src/HouseHolder/dense.h"

TEST(QR, B) {

    std::vector<double> a = { 3, 2, 4, 4};
    int n = sqrt(a.size());
    Matr Q(n), R(n);
    Matr A(a);
    A.out();
    Q = A.QRH().first;
    Matr B(a);
    R = B.QRH().second;
    Q.out();
    R.out();
    //Householder method

    Matr C(a);
    Q = C.QRG().first;
    Q.out();
    Matr D(a);
    R = D.QRG().second;
    R.out();

}
TEST(QR, T) {
    std::vector<double> a = { 3, 1,4,4, 5, 3, 1, 2, 3};
    Matr A(a);
    A.out();
    A.QRH().first.out();
    Matr B(a);
    B.QRH().second.out();
//householder working

    Matr C(a);
    C.out();
    C.QRG().first.out();
    Matr D(a);
    D.QRG().second.out();

    //All is working;)
}

TEST (GMRES, ZHEST) {
    std::vector<double> a = {3, 1, 4, 2, 5, 3, 3, 1, 3};
    std::vector<double> x0 = {1,1,1};
    Matr A(a);
    A.out();
    std::vector<double> b = {1, 2, 3};
    std::vector<double> x = A.GMRES(A,x0,b,0.00000001,3);
    for (int i = 0; i < 3; i++) {
        std::cout << x[i] << "\n";
    }
}
TEST (GMRES, FULLZHEST) {
    std::vector<double> a = {3, 1, 0, 2, 1, 3, 0, 1, 3, 0,5,6,7,1,1,0};
    std::vector<double> x0 = {1,1,1,1};
    Matr A(a);
    A.out();
    std::vector<double> b = {1, 2, 3,4};
    std::vector<double> x = A.GMRES(A,x0,b,0.00000001,3);
    for (int i = 0; i < 4; i++) {
        std::cout << x[i] << "\n";
    }
}

TEST (BCG, 1) {
    std::vector<double> a = {3, 1, 0, 2, 1, 3, 0, 1, 3, 0,5,6,7,1,1,0};
    std::vector<double> x0 = {1,1,1,1};
    Matr A(a);
    A.out();
    std::vector<double> b = {1, 2, 3,4};
    std::vector<double> x = A.GMRES(A,x0,b,0.00000001,3);
    for (int i = 0; i < 4; i++) {
        std::cout << x[i] << "\n";
    }
}

TEST (CGS, 1) {
    std::vector<double> m =  {1,2,3,3,2,-1,2,1,0.9}, x0 = {1,1,1}, b = {1,2,3};
    Matr M(m);
    M.out();
    std::vector<double> result = M.CGS(M,x0,b,0.000000001);
    for (int i = 0; i < x0.size(); i++) {
        std::cout << result[i] << "\n";
    }
    //god damn, hard one
}