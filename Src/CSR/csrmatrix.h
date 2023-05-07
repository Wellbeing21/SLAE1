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


        Csr_matrix(int n);

        std::vector<double> data;
        std::vector<std::size_t> col_ind;
        std::vector<std::size_t> sup_row;


    public:

        double operator()(std::size_t i, std::size_t j) const;
        void out(Csr_matrix &A);
        void Wolout(Csr_matrix &A);


        std::vector<double> operator*(const std::vector<double> &fre) const;
        std::vector<double> ReshIT(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b, const double &t, const double &e);
        std::vector<double> ReshITS(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,const int &rd,const  double &e,const double &lam1,const double &lam2);
        std::vector<double> ReshGZ(const Csr_matrix &A, const std::vector<double> &x01,const std::vector<double> &b,const  double &e);
        std::vector<double> ReshYak(const Csr_matrix &A,const  std::vector<double> &x01, const std::vector<double> &b,const  double &e);
        std::vector<double> ReshGZS(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b,const double &e);
        std::vector<double> ReshSOR(const Csr_matrix &A, const std::vector<double> &x01,const std::vector<double> &b,const double &e,const double &w);
        std::vector<double> ReshSSOR(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,const  double &e, const double &w);
        std::vector<double> ReshSSORCHEB(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b, const double &e, const double &p, const double &w);
        double eigenvalue(const Csr_matrix &M, const std::vector<double> &x01,const double &e);
        std::vector<double> ReshRapidGradient(const Csr_matrix &A,const std::vector<double> &x01,const std::vector<double> &b, const double &e);
        std::vector<double> ReshNesterov(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,const double &e);
        double Amult(const Csr_matrix &A,const std::vector<double> &x2);
        std::vector<double> ReshConjugateGradient(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b,const double &e);
        Csr_matrix(std::vector<DOK> vec);
        Csr_matrix(std::vector<double> vec);

        void LU0(const Csr_matrix &A);
        void LU(const Csr_matrix &A);

        Csr_matrix CHOL0(const Csr_matrix &A);

        std::vector<double>
        ReshConjugateGradientObusl(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,
                                   const double &e);

        std::vector<double> SolverGaussReversed_forDownTriang(const Csr_matrix &A, const std::vector<double> &b);
    };

};




#endif