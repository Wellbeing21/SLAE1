#ifndef SLAE_CSRMATRIX_H_GLEB
#define SLAE_CSRMATRIX_H_GLEB
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
struct DOK {
    int i;
    int j;
    double val;
};
double mistakemod(const std::vector<double>&x1, const std::vector<double>&x2);
double vecmod(std::vector<double>&x1);
double scalarmult(std::vector<double>&x1, std::vector<double>&x2);

namespace CSR_space {
    class Csr_matrix {
    private:


    public:
        std::vector<double> data;
        std::vector<std::size_t> col_ind;
        std::vector<std::size_t> sup_row;

        double operator()(std::size_t i, std::size_t j);
        void out(Csr_matrix &A);
        void Wolout(Csr_matrix &A);


        std::vector<double> operator*(const std::vector<double> fre);
        std::vector<double> ReshIT(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double t, double e);
        std::vector<double> ReshITS(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, int rd, double e, double lam1, double lam2);
        std::vector<double> ReshGZ(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        std::vector<double> ReshYak(Csr_matrix &A,  std::vector<double> x0, const std::vector<double> b, double e);
        std::vector<double> ReshGZS(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        std::vector<double> ReshSOR(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e,double w);
        std::vector<double> ReshSSOR(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e,double w);
        std::vector<double> ReshSSORCHEB(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e,double p, double w);
        double eigenvalue(Csr_matrix &M, std::vector<double> &x0, double e);
        std::vector<double> ReshRapidGradient(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        std::vector<double> ReshNesterov(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        double Amult(Csr_matrix &A, std::vector<double>&x2);
        std::vector<double> ReshConjugateGradient(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        Csr_matrix(std::vector<DOK> vec);

    };

};




#endif