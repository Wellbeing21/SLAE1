#ifndef SLAE_CSRMATRIX_H_GLEB
#define SLAE_CSRMATRIX_H_GLEB
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>




class Matr {
private:
    std::vector <double> vv;
public:
    Matr(std::vector<std::vector<double>> &a);
    Matr(int n);
    Matr (std::vector <double> &a);
    void out();
    double operator()(int i, int j) const;
//with const if we don't need to change something
    double operator()(int i, int j) ;
//without(to change smth)
    Matr operator* (Matr &M);
    std::pair<Matr, Matr> QRH ();
    std::pair<Matr, Matr> QRG ();

    std::vector<double>  operator* (const std::vector<double> &a) const;
    std::vector<double> SolverGaussReversed_forUpTriang(const Matr A, const std::vector<double> b);
    std::vector<double> GMRES ( const Matr &A,  const std::vector<double> &x01, const std::vector<double> &b,const double e,const int m1);
    std::vector<double> TransposedMatrixMult (const std::vector<double> &a) const;
    std::vector<double> CGS ( const Matr &A,  const std::vector<double> &x0, const std::vector<double> &b,const double e);
    std::vector<double> BCG ( const Matr &A,  const std::vector<double> &x0, const std::vector<double> &b,const double e);

};


#endif