#include <gtest/gtest.h>
#include "../Src/solvers/solve_tridiagonal.hpp"
#include "../Src/CSR/csrmatrix.h"
using namespace CSR_space;
TEST(A, B){
    std::vector<line> Mat = {{0, 2, 1},{1,10,-5},{1,-5,2},{1,4,0}};
    std::vector<double> free1 = {-5, -18, -40, -27};
    TriMat matr1;
    matr1.apply(Mat);
    std::vector<double> solut = solver::findsolution(matr1, free1);
    for (auto i: solut) {
        std::cout << i << " ";
    }
}
TEST (B, C){
    std::vector<double> b;
    std::vector<DOK> a;
    a.resize(4);
    //there's massive: 1    5
    //                 0    2
    //                 0    0
    //                 0    1
    a[0] = {0, 0, 1};
    a[1] = {0, 1, 5};
    a[2] ={1, 1, 2};
    a[3] = {3, 1, 1};
    CSR_space::Csr_matrix matrix(a);
    std::cout << matrix(3,1) << "\n";
    //1-stt element
    std::cout << matrix(1,1) << "\n";
    //second
    std::cout << matrix(2,0) << "\n";
    //third
    b.resize(2);
    b = {1, 5};
    b = matrix*b;
    for (auto el : b) {
        std::cout << el << " ";
    }

}
TEST (YAK, SOS){
    double e = 0.000000000000001;
    //test for yakobi, its working
    std::cout << "Test for Yakobi:\n";
    std::vector<DOK> AYAK {{0,0,0.7},{0,1,0.2}, {1,1,0.2}, {2,0,0.1}, {2,2,0.3} };
    Csr_matrix MYAK(AYAK);
    MYAK.out(MYAK);
    std::vector<double> x0YAK{1, 1, 1}, bYAK{1,1,1};
    //taking vector, now we need to use yakobi method
    x0YAK = MYAK.ReshYak(MYAK,x0YAK,bYAK,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x0YAK[i] << "\n";
    }

}
TEST (GZ, SOS){
    double e = 0.000001;
    //it's test for my G-Z method
    std::cout << "Test for Gauss-Zeidel:\n";
    std::vector<DOK> AGZ {{0,0,4},{0,1,1},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix MGZ(AGZ);
    MGZ.out(MGZ);
    std::vector<double> x0GZ{1, 1, 3}, bGZ{1,1,1};
    x0GZ = MGZ.ReshGZ(MGZ,x0GZ,bGZ,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x0GZ[i] << "\n";
    }

}
TEST (IT, SOS){
    double e = 0.000000000001;
    double t = 0.1;
    //it's test for my iterations method
    std::cout << "Test for Simple iterations(samuj sosovuj):\n";
    std::vector<DOK> Ait {{0,0,3},{0,1,0.5},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix Mit(Ait);
    Mit.out(Mit);
    std::vector<double> x0it{1, 1, 1}, bit{1,1,1};
    x0it= Mit.ReshIT(Mit,x0it,bit,t,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x0it[i] << "\n";
    }
}
TEST (IT, NESOS){
    double e = 0.00000000000000000000000000000000000001;
//Test for Smart iterations (with polinoms)
    std::cout << "Test for Smart iterations (with polinoms):\n";
    std::vector<DOK> AITS {{0,0,3},{0,1,0.5},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix MITS(AITS);
    MITS.out(MITS);
    std::vector<double> x0ITS{1, 1, 1}, bITS{1,1,1};
    x0ITS= MITS.ReshITS(MITS,x0ITS,bITS,3,e,2,8);
    for (int i = 0; i < 3; i++) {
        std::cout << x0ITS[i] << "\n";
    }
    std::cout << "Really ne sos, see that accuracy;)\n";


}

TEST(GZS, SOVSEMSOS) {
    double e = 0.00000001;
    std::cout << "Test for GZ(symmetrical):\n";
    std::vector<DOK> AGZS {{0,0,3},{0,1,1},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix MGZS(AGZS);
    MGZS.out(MGZS);
    std::vector<double> x0GZS{1, 1, 3}, bGZS{1,1,1};
    x0GZS = MGZS.ReshGZS(MGZS,x0GZS,bGZS,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x0GZS[i] << "\n";
    }
    //tested Zeydel Symmetrical

}
TEST(SOR, SOS) {
    double e = 0.0000000001;
    double w, mu = (double)2 / 3;// mu = max eigenvalue from C = I - D^-1 A
    w = 1 + (mu * mu) / (1 + sqrt(1 - mu * mu));
    //tests for SOR:
    std::cout <<" Test for SOR:\n";
    std::vector<DOK> SORv {{0,0,1},{0,1,-2},{1,1,1},{1,2,1},{2,1,1}, {2,2,3} };
    Csr_matrix MSOR(SORv);
    MSOR.out(MSOR);
    std::vector<double> xsor{1, 1, 3}, bsor{1,1,1};
    xsor = MSOR.ReshSOR(MSOR,xsor,bsor,e,w);
    for (int i = 0; i < 3; i++) {
        std::cout << xsor[i] << "\n";
    }
}

TEST(SSOR, SSOS) {
    double e = 0.000000000001;
    double w, mu = (double)2 / 3;// mu = max eigenvalue from C = I - D^-1 A
    w = 1 + (mu * mu) / (1 + sqrt(1 - mu * mu));
    //tested Zeydel Symmetrical
    //Tests for SSOR:
    std::cout << "Test for SSOR:\n";
    std::vector<DOK> SSORv {{0,0,1},{0,1,-2},{1,1,1},{1,2,1},{2,1,1}, {2,2,3} };
    Csr_matrix MSSOR(SSORv);
    MSSOR.out(MSSOR);
    std::vector<double> xssor{1, 1, 3}, bssor{1,1,1};
    xssor = MSSOR.ReshSSOR(MSSOR,xssor,bssor,e,w);
    for (int i = 0; i < 3; i++) {
        std::cout << xssor[i] << "\n";
    }

}

TEST(SSORCHEB, SOVSEMNESOS) {
    double e = 0.00000001;
    double w, mu = (double)2 / 3;// mu = max eigenvalue from C = I - D^-1 A
    w = 1 + (mu * mu) / (1 + sqrt(1 - mu * mu));
    //Tests for SSOR with cheb:
    std::cout << "Test for SSOR with cheb. shit:\n";
    std::vector<DOK> SSORchv {{0,0,1},{0,1,-2},{1,1,1},{1,2,1},{2,1,1}, {2,2,3} };
    Csr_matrix MSSORch(SSORchv);
    MSSORch.out(MSSORch);
    std::vector<double> xssorch{1, 1, 3}, bssorch{1,1,1};
    xssorch = MSSORch.ReshSSORCHEB(MSSORch,xssorch,bssorch,e,mu, w);
    for (int i = 0; i < 3; i++) {
        std::cout << xssorch[i] << "\n";
    }
}

TEST(eigenvector, eigenvalue) {
    //Tests for max eigenvalue
    double e = 0.0000000001;
    std::cout << "Test for max value:\n";
    std::vector<DOK> maxv {{0,0,1},{0,1,-2},{1,1,1},{1,2,1},{2,1,1}, {2,2,5} };
    Csr_matrix maxM(maxv);
    maxM.out(maxM);
    std::vector<double> xmax{1, 1, 3};
    double a = maxM.eigenvalue(maxM,xmax,e);
    std::cout << a << " - max eigenvalue\n";
    for (int i = 0; i < 3; i++) {
        std::cout << xmax[i] << "\n";
    }
}

TEST (GRAD, SOS) {
    //Tests for Gradient rapid thing
    double e = 0.000000001;
    std::cout << "Test for gradient test(naiskoreyshij):\n";
    std::vector<DOK> grrapv {{0,0,0.4},{1,1,6},{1,2,0.1},{2,1,0.1}, {2,2,0.5} };
    Csr_matrix grrapidM(grrapv);
    grrapidM.out(grrapidM);
    std::vector<double> xgrrapid{2, 1, 3}, bgrrapid{2,5,1};
    xgrrapid = grrapidM.ReshRapidGradient(grrapidM,xgrrapid,bgrrapid,e);
    for (int i = 0; i < 3; i++) {
        std::cout << xgrrapid[i] << "\n";
    }


}

TEST (GRADNEST, NESOS) {
    //Tests for Gradient NESTEROV
    double e = 0.000000001;
    std::cout << "Test for gradient Nesterov:\n";
    std::vector<DOK> nesv {{0,0,0.4},{1,1,6},{1,2,0.1},{2,1,0.1}, {2,2,0.5} };
    Csr_matrix nesM(nesv);
    nesM.out(nesM);
    std::vector<double> xnes{1, 1, 1}, bnes{2,5,1};
    xnes = nesM.ReshNesterov(nesM,xnes,bnes,e);
    for (int i = 0; i < 3; i++) {
        std::cout << xnes[i] << "\n";
    }
}
TEST (GRADCONJUGATE, ABSNESOS) {
    //Tests for Gradient NESTEROV
    double e = 0.000000000001;
    std::cout << "Test for Conjugate gradient method(ABSOLUTELY NE SOS):\n";
    std::vector<DOK> conjv {{0,0,0.5},{1,1,0.4},{1,2,0.3},{2,1,0.3}, {2,2,0.7} , {3,3,1}, {4,4,0.4}, {4,5,3},{5,5,3.2},{6,6,11},{7,6,0.1},{7,7,{0.6}},{8,2, 0.3},{8,8,1.3},{9,9,1}};
    Csr_matrix conjM(conjv);
    conjM.out(conjM);
    std::vector<double> xconj{1, 1, 1,1,1,1,1,1,1,1}, bconj{2,5,1,3,1,1,1,1,1,1};
    xconj = conjM.ReshConjugateGradient(conjM,xconj,bconj,e);
    for (int i = 0; i < xconj.size(); i++) {
        std::cout << xconj[i] << "\n";
    }
}

TEST (KR2,1) {
    int n = 289;
    double e = pow(10,-13);
    double t = 0.1;
    double b = 40, a = 18, c = 2;
    //it's test for my iterations method
    std::cout << "Test for Simple iterations(samuj sosovuj):\n";
    std::vector<DOK> A(n+(n-1)*2 + (n- sqrt(n))*2);
    int it = 0;
    std::vector<double> c1(n);

    A[it] = {0,0, 2*b};
    it++;
    A[it] = {0, 1, a};
    it++;
    A[it] = {0, sqrt(n), a};
    it++;
    //First is ready
    for (int i = 1; i < n-1; i++) {
        if (i > sqrt(n)-1) {
            A[it] = {i, i - sqrt(n), a};
            it++;
        }
        A[it] = {i,i - 1, a};
        it++;
        A[it] = {i,i,2*b};
        it++;
        A[it] = {i, i+1, a};
        it++;
        if (i < n - sqrt(n)) {
            A[it] = {i, i+ sqrt(n), a};
            it++;
        }
    }
    //Fill n-1 strings
    A[it] = {n-1,n-1- sqrt(n), a };
    it++;
    A[it] = {n-1,n-2, a};
    it++;
    A[it] = {n-1,n-1, 2*b};
    it++;
    Csr_matrix M(A);
    //M.out(M);
    //its working
    for (int i = 0; i < n; i++){
        c1[i] = c;
    }
    double lammax = 2 * (b + 2*a*cos(atan(1)*4/(sqrt(n)+1)));
    double lammin = 2 * (b - 2*a*cos(atan(1)*4/(sqrt(n)+1)));
    std::vector<double> x0(289);
    x0 = M.ReshIT(M,x0,c1,1/lammax,e);
    std::cout << "That was 1)a)\n";
    for (int i = 0; i<n; i++) {
        x0[i] = 0;
    }

    x0 = M.ReshRapidGradient(M,x0,c1,e);
    for (int i =0; i<n; i++) {
        x0[i] =0;
    }
    std::cout << "That was 1)b)Rapid gradient\n";
    x0 = M.ReshITS(M,x0,c1,3,e,lammin*0.9, lammax*1.1);
    for (int i =0; i<n; i++) {
        x0[i] =0;
    }
    std::cout << "That was 1)c)Chebushov\n";
    //we need to know maximum eiegenvalue of C, so i will check it
    it = 0;

    double w = 0.5;
    x0 = M.ReshSSOR(M,x0,c1, e,w);
    for (int i = 0; i<n; i++) {
        x0[i] =0;
    }
    std::cout << "That was 1)d)Symmetrical \n";
}

TEST (KR2,2) {
    int n = 4;
    double e = pow(10,-13);
    std::vector<double> x0(4);
    std::vector<DOK> A {{0,0,12},{1,1,14},{2,2,16},{3,3,18}};
    Csr_matrix M(A);
    double lammin = 12, lammax = 18;
    //eigenmin = (1,0,0,0); eigenmax = (0,0,0,1);
    std::vector<double> b(4);
    for (int i = 0; i < n; i++) {
        b[i] = 4;
    }
    double tau = 0.9 * 2 / lammax;
    x0 = M.ReshIT(M,x0,b,tau,e);
    std::cout << "That was 2)a)with tau = 0.9*2 / lammax \n";
    for (int i = 0; i < n; i++) {
        x0[i] = 0;
    }
    x0 = M.ReshRapidGradient(M,x0,b,e);
    std::cout << "That was 2)b)with opt tau \n";
    for (int i = 0; i < n; i++) {
        x0[i] = 0;
    }
    x0 = M.ReshITS(M,x0,b,3,e,lammin,lammax);
    std::cout << "That was 2)c)with opt tau \n";
    for (int i = 0; i < n; i++) {
        x0[i] = 0;
    }
    x0 = M.ReshConjugateGradient(M,x0,b,e);



}
TEST (test_LU, 1) {
    std::vector<double> a = {5,1,0,0,1,  0,5,1,0,1, 1,1,5,1,0,  0,0,1,5,0,   0,1,1,1,5};
    Csr_matrix M(a);
    M.out(M);
    M.LU(M);
    std::cout << "That's the LU-decomposition\n";
    M.LU0(M);
    std::cout << "That's the LU-decomposition(zero)\n";
    std::cout << "You can add it to Wolfram Alpha just ctrl + c ...";
}

TEST (test_CHOL, 1) {
    std::vector<double> a = {6,1,0,0,1,  1,5,1,0,1, 0,1,5,1,0,  0,0,1,3,0,   1,1,0,0,5};
    Csr_matrix M(a);
    M.Wolout(M);
    M.out(M);
    M.CHOL0(M);
}
TEST (test_CHOL, 2) {
    std::vector<double> a = {5,0,1,  0,4,1,  1,1,7};
    Csr_matrix M(a);
    M.Wolout(M);
    M.out(M);
    Csr_matrix A = M.CHOL0(M);
    std::cout << "You can check it in Wolfram\n";
}

TEST (Gauss_triag, 1) {
    std::vector<double> a = {5,0,0,  3,4,0,  1,1,7};
    Csr_matrix M(a);
    M.Wolout(M);
    std::vector<double> b = {1,1,1};
    std::vector<double> x(3);
    M.out(M);
    x = M.SolverGaussReversed_forDownTriang(M, b);
    for (int i = 0; i < 3; i++){
        std::cout << x[i] <<"\n";
    }

}


TEST (CG_obusl, 1) {
    double e = 0.003;
    int n = 5;
    std::vector<double> a = {6,1,0,0,1,  1,5,1,0,1, 0,1,5,1,0,  0,0,1,3,0,   1,1,0,0,5};
    Csr_matrix M(a);

    std::vector<double> b = {1,1,1,1,1}, x0 = {1,1,1,1,1};
    M.Wolout(M);
    M.out(M);
    std::vector<double> x(n);
    x = M.ReshConjugateGradientObusl(M,x0,b,e);
    for (int i = 0; i < n; i++){
        std::cout << x[i] << "\n";
    }

}

TEST (CG_obusl, 2) {
    double e = 0.0001;
    int n = 3;
    std::vector<double> a = {6,1,0,  1,5,1, 0,1,5};
    Csr_matrix M(a);

    std::vector<double> b = {1,1,1}, x0 = {1,1,1};
    M.Wolout(M);
    M.out(M);
    std::vector<double> x(n);
    x = M.ReshConjugateGradientObusl(M,x0,b,e);
    for (int i = 0; i < n; i++){
        std::cout << x[i] << "\n";
    }

}