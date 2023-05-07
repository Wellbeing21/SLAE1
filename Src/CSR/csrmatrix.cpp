//
// Created by gorilla on 18.02.23.
//

#include "csrmatrix.h"

double mistakemod(const std::vector<double>&x1,const std::vector<double>&x2) {
    double a = 0;
    for (int i = 0; i < x1.size(); i++) {
        a += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
      a = std::sqrt(a);
    return a;
}

double vecmod(std::vector<double>&x1) {
    double a = 0;
    for (int i = 0; i < x1.size(); i++) {
        a += (x1[i]) * (x1[i]);
    }
    a = std::sqrt(a);
    return a;
}

double operator* (const std::vector<double> &a, const std::vector<double> &b) {
    double sum = 0;
    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<double> operator* (const double &a, const std::vector<double> &b) {
    std::vector<double> res(b.size());
    for (int i = 0; i < b.size(); i++) {
        res[i] = a * b[i];
    }
    return res;
}

std::vector<double> operator+ (const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(b.size());
    for (int i = 0; i < b.size(); i++) {
        res[i] = a[i] + b[i];
    }
    return res;
}

std::vector<double> operator- (const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> res(b.size());
    for (int i = 0; i < b.size(); i++) {
        res[i] = a[i] - b[i];
    }
    return res;
}


double CSR_space::Csr_matrix::operator()(std::size_t i, std::size_t j) const{
    int c = 0;
    for (std::size_t it = sup_row[i]; it < sup_row[i+1]; it++) {
        if (col_ind[it] == j) {
            return data[it];
            c = 1;
        }}
        if (c == 0) {
            return 0;
    }
}

std::vector<double> CSR_space::Csr_matrix::operator*( const std::vector<double> &fre) const{
    std::vector<double> a;
    std::size_t n = sup_row.size() - 1;
    a.resize(n);
    for (int it = 0; it < n; it++) {
        for (int jt = sup_row[it]; jt < sup_row[it+1]; jt++) {
            a[it] += data[jt] * fre[col_ind[jt]];
        }
    }
    return a;
}

void CSR_space::Csr_matrix::out(Csr_matrix &A) {
    int n = (A.sup_row.size()-1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout  << std::setw(6) << A(i,j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}
void CSR_space::Csr_matrix::Wolout(Csr_matrix &A){
    int n = A.sup_row.size()-1;
    std::cout << "{";
    for (int i = 0; i < n; i++) {
        std::cout << "{";
        for (int j = 0; j < n; j++) {
            if (j != n - 1) {std::cout<<A(i,j) << ", ";}
            else {std::cout<<A(i,j);}
        }
        if (i != n - 1) {std::cout << "}, ";}
        else {std::cout << "}";}
    }
    std::cout << "}\n";
}

std::vector<double> CSR_space::Csr_matrix::ReshYak(const Csr_matrix &A,const  std::vector<double> &x01, const std::vector<double> &b,const  double &e) {
    std::vector<double> x0 = x01;
    int n = (A.sup_row.size() - 1);
    std::vector<double> ans;
    double r = 1;
    //multiplied without diagonal
    while (not(r < e)) {
        double d;
        for (int i = 0; i < n; i++) {
            x0[i] = b[i];
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {
                    d = data[s];
                }
                else {
                    x0[i] -= x0[col_ind[s]] * data[s];
                }
            }
            x0[i] /= d;
        }
        ans = A * x0;
        r = mistakemod(ans, b);
    }
    //count our neuvyazky po supremumu
    return x0;

}

std::vector<double> CSR_space::Csr_matrix::ReshIT(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,const double &t, const double &e) {
    int n = x01.size();
    std::vector<double> x0 = x01;
    int iter = 0;
    std::vector<double> mul(n), bf(n);
    double r = 1;
    while (std::fabs(r) > e) {
        mul = A * x0;
        iter++;
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] - t * mul[i] + t * b[i];
        }
        bf = A * x0;
        r = mistakemod(bf, b);
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;


}

//same RESH with simple iterations, but with smart polinoms
//n = degree of our polynom
std::vector<double> CSR_space::Csr_matrix::ReshITS(const Csr_matrix &A,const std::vector<double> &x01,const  std::vector<double> &b, const int &rd,const double &e,const double &lam1,const double &lam2) {
    int n = A.sup_row.size() - 1;
    std::vector<double> x0 = x01;
    int num = std::pow(2, rd);
    std::vector<double> mul, bf, t(num), ord(num), x(num);
    //first, we need to make our t-vector with cycle with num roots
    double sina = sin(4 * atan(1) / num), cosa = cos(4 * atan(1) / num), sinb = sin(4 * atan(1) / (2 * num));
    x[0] = cos(4 * atan(1) / (2 * num));
    t[0] = 1 / ((lam1 + lam2) / 2 + x[0] * (lam2 - lam1) / 2);
    for (int i = 1; i < num; i++) {
        x[i] = x[i - 1] * cosa - sinb * sina;
        t[i] = 1 / ((lam1 + lam2) / 2 + x[i] * (lam2 - lam1) / 2);
        sinb = sinb * cosa + x[i - 1] * sina;
    }
    //ok we found xk and tk
    //now we need order here to take tk with it's number
    double j = num / 2;
    ord[0] = 0;
    ord[j] = 1;
    for (int i = 1; i < n; i++) {
        j = j/2;
        if (j < 1) {break;}
        for (int jt = j; jt < num; jt = jt + 2*j) {
            ord[jt] = std::pow(2,i+1) - ord[jt - j] - 1;
        }
    }
    //and we found our order so, let's start iterations
    double r = 1;
    int por = 0;
    int iter = 0;
    while (std::fabs(r) > e) {
        iter++;
        mul = A*x0;
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] - t[ord[por]] * mul[i] + t[ord[por]] * b[i];
        }
        por++;
        if (por >= num-1) {
            por = 0;
        }
        //we increase our t_por
        bf = A*x0;
        r = mistakemod(bf, b);


    }
    std::cout << iter <<"- this is the number of iterations" << "\n";
    return x0;

}


///i changed Zeidel here
std::vector<double> CSR_space::Csr_matrix::ReshGZ(const Csr_matrix &A, const std::vector<double> &x01,const std::vector<double> &b,const  double &e) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() -1;
    double r = 1;
    std::vector<double> bf;
    while (not r < e) {
        double d = 0;
        for (int i = 0; i < n; i++) {
            x0[i] = b[i];
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {d = data[s];}
                else {
                    x0[i] -= x0[col_ind[s]] * data[s];
                }
            }
            x0[i] /= d;
        }
        bf = A*x0;

        r = mistakemod(bf,b);
        //finded supremum
    }
    return x0;
}
//Gauss zeydel symmetrical double time (x_{i+05}, x_{i+1}}
std::vector<double> CSR_space::Csr_matrix::ReshGZS(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b,const double &e) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() -1;
    double r = 1;
    std::vector<double> bf;
    while (not r < e) {
        double d = 0;
        for (int i = 0; i < n; i++) {
            x0[i] = b[i];
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {d = data[s];}
                else {
                    x0[i] -= x0[col_ind[s]] * data[s];
                }
            }
            x0[i] /= d;
        }
        //xi+0.5 done
        for (int i = n-1; i >= 0; i--) {
            x0[i] = b[i];
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {d = data[s];}
                else {
                    x0[i] -= x0[col_ind[s]] * data[s];
                }
            }
            x0[i] /= d;
        }
        bf = A*x0;
        r = mistakemod(bf,b);
        //finded supremum
    }
    return x0;
}
///sor method upper relaxation
std::vector<double> CSR_space::Csr_matrix::ReshSOR(const Csr_matrix &A, const std::vector<double> &x01,const std::vector<double> &b,const double &e,const double &w) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() -1;
    double r = 1;
    std::vector<double> bf;
    while (not r < e) {
        double d = 0, temp = 0;
        for (int i = 0; i < n; i++) {
            temp = w*b[i]; // because we can't change x0[i] - we'll multiply it on (w-1)D
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {
                    d = data[s];
                    temp -= (w-1) * x0[col_ind[s]] * data[s];
                }
                else {
                    temp -= w * x0[col_ind[s]] * data[s];
                }
            }
            temp /= d;
            x0[i] = temp;
        }
        bf = A*x0;

        r = mistakemod(bf,b);
        //finded supremum
    }
    return x0;
}
///plus one
std::vector<double> CSR_space::Csr_matrix::ReshSSOR(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,const  double &e, const double &w) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() -1;
    double r = 1;
    std::vector<double> bf;
    int iter = 0;
    while (not (r < e)) {
        iter++;
        double d = 0, temp = 0;
        for (int i = 0; i < n; i++) {
            temp = w*b[i]; // because we can't change x0[i] - we'll multiply it on (w-1)D
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {
                    d = data[s];
                    temp -= (w-1) * x0[col_ind[s]] * data[s];
                }
                else {
                    temp -= w * x0[col_ind[s]] * data[s];
                }
            }
            temp /= d;
            x0[i] = temp;
        }
        //xi+0.5 done
        for (int i = n-1; i >= 0; i--) {
            temp = w*b[i];
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {
                    d = data[s];
                    temp -= (w-1) * x0[col_ind[s]] * data[s];
                }
                else {
                    temp -= w * x0[col_ind[s]] * data[s];
                }
            }
            temp /= d;
            x0[i] = temp;
        }
        bf = A*x0;

        r = mistakemod(bf,b);



        //finded supremum
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;
}

///And now chebushevskoe for SSOR
std::vector<double> CSR_space::Csr_matrix::ReshSSORCHEB(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b, const double &e, const double &p, const double &w) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() -1;
    double r = 1;
    std::vector<double> bf, mu(3), ypred = x0;
    mu[0] = 1;
    mu[1] = 1/p;
    //p - spectral radius
    //
    double iter = 0;
    while (not (r < e)) {
        iter++;
        mu[2] = 2 / p * mu[1] - mu[0];
        double c1 = 2 * mu[1] / (p * mu[2]), c2 = mu[0] / mu[2];
        double d = 0, temp = 0;
        for (int i = 0; i < n; i++) {
            temp = w*b[i]; // because we can't change x0[i] - we'll multiply it on (w-1)D
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {
                    d = data[s];
                    temp -= (w-1) * x0[col_ind[s]] * data[s];
                }
                else {
                    temp -= w * x0[col_ind[s]] * data[s];
                }
            }
            temp /= d;
            x0[i] = temp;
        }
        //xi+0.5 done
        for (int i = n-1; i >= 0; i--) {
            temp = w*b[i];
            for (int s = sup_row[i]; s < sup_row[i + 1]; s++) {
                if (col_ind[s] == i) {
                    d = data[s];
                    temp -= (w-1) * x0[col_ind[s]] * data[s];
                }
                else {
                    temp -= w * x0[col_ind[s]] * data[s];
                }
            }
            temp /= d;
            x0[i] = temp;

        }
        //we've finded Py+c now c1 * ... + c2 * ypred; now yi+1
        for (int i = 0; i < n; i++) {
            x0[i] = c1 * x0[i] - c2 * ypred[i];
        }
        mu[0] = mu[1];
        mu[1] = mu[2];
        ypred = x0;
        bf = A*x0;

        r = mistakemod(bf,b);
        //finded supremum
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;
}

double CSR_space::Csr_matrix::eigenvalue(const Csr_matrix &M, const std::vector<double> &x01,const double &e) {
    std::vector<double> x0 = x01;
    for (int i = 0; i < x0.size(); i++) {
        x0[i] = x0[i] / vecmod(x0);
    }
    std::vector<double> xpred;
    double r = 1;
    while (not(r < e/10000)) {
        xpred = x0;
        x0 = M*x0;
        for (int i = 0; i < x0.size(); i++) {
            x0[i] = x0[i] / vecmod(x0);
        }
        r = mistakemod(x0, xpred);
    }
    xpred = M*x0;
    return vecmod(xpred);
}

std::vector<double> CSR_space::Csr_matrix::ReshRapidGradient(const Csr_matrix &A,const std::vector<double> &x01,const std::vector<double> &b, const double &e) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() - 1;
    int iter = 0;
    std::vector<double> mul, bf, ri(n);
    double r = 1, a;
    while (std::fabs(r) < e) {
        mul = A * x0;
        for (int i = 0; i < n; i++ ) {
            ri[i] = mul[i] - b[i];  //finded nevyazku
        }
        std::vector<double> znam;
        znam = A*ri;
        a = (ri*ri) / (ri*znam);
        //got alpha, now do iterations
        iter++;
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] - a * ri[i];
        }
        r = vecmod(ri);
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;

}
double CSR_space::Csr_matrix::Amult(const Csr_matrix &A,const std::vector<double> &x2) {
    std::vector<double> x0 = x2;
    double sc =0;
    std::vector x11 = x2;
    std::vector x22 = x2;
    x22 = A*x2;
    for (int i = 0; i < x11.size(); i++) {
        sc += x11[i] * x22[i];
    }
    return sc;
}

std::vector<double> CSR_space::Csr_matrix::ReshNesterov(const Csr_matrix &A, const std::vector<double> &x01, const std::vector<double> &b,const double &e) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() - 1;
    int iter = 0;
    std::vector<double> mul, bf, ri(n), xpr = x0 , xnow(n);
    double r = 1, a;
    while ( std::fabs(r) > e) {
        mul = A * x0;
        for (int i = 0; i < n; i++ ) {
            ri[i] = mul[i] - b[i];  //finded nevyazku
        }
        std::vector<double> znam;
        znam = A*ri;
        a = (ri*ri) / (ri*znam);
        //got alpha, now do iterations
        iter++;
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] - a * ri[i];
        }
        xnow = x0;
        //do gradient thing
        //want to find y
        for (int i =0; i < n; i++) {
            x0[i] = x0[i] + iter / (iter + 3) * (x0[i] - xpr[i]);
        }
        xpr = xnow;
        r = vecmod(ri);
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;

}
std::vector<double> CSR_space::Csr_matrix::ReshConjugateGradient(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b,const double &e) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() - 1;
    int iter = 0;
    std::vector<double> Ax, di(n), di1(n), ri(n), ri1(n);
    double r = 1, alpha;
    //need to do 1-st iteration
    Ax = A*x0;
    ri1 = Ax - b;
    di1 = ri1;
    alpha = (di1*ri1) / Amult(A,di1);
    //got alpha, now do iteration
    iter++;
    x0 = x0 - (alpha * di1);
    ri = ri1;
    di = di1;
    r = vecmod(ri);
    while (true) {
        Ax = A * x0;
        ri1 = Ax - b;  //finded nevyazku(r_{i+1})
        if (vecmod(ri1) < e) {break;}
        double beta = (ri1*ri1) / (di*ri);
        di1 = ri1 + (beta * di);
        //finded d_{i+1}
        alpha = (di1*ri1) / Amult(A,di1);
        //got alpha, now do iteration
        iter++;
            x0 = x0 - (alpha * di1);
        ri = ri1;
        di = di1;
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;

}


void CSR_space::Csr_matrix::LU(const Csr_matrix &A) {
    int n = A.sup_row.size()-1;
    std::vector<double> l(n * n);
    std::vector<double> u(n * n);
    for (int i = 0; i < n; i ++) {
        l[i * n + i] = 1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) {
                u[i * n + j] = A(i,j);
                for (int k = 0; k < i; k++) {
                    u[i * n + j] -= l[i * n + k] * u[k * n + j];
                }
            }
            else {
                l[i * n + j] = A(i,j);
                for (int k = 0; k < j; k++) {
                    l[i * n + j] -= l[i * n + k] * u[k * n + j];
                }
                l[i * n + j] = l[i * n + j]/u[j * n + j];
            }
        }

    }

    Csr_matrix L(l);
    Csr_matrix U(u);
    L.Wolout(L);
    U.Wolout(U);
}

void CSR_space::Csr_matrix::LU0(const Csr_matrix &A) {
    int n = A.sup_row.size()-1;
    std::vector<double> l(n * n);
    std::vector<double> u(n * n);
    std::vector<double> lvl(n * n);
    for (int i = 0; i < n; i ++) {
        l[i * n + i] = 1;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A(i,j) != 0) {
                lvl[i * n + j] = 1;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) {
                u[i * n + j] = A(i,j);
                for (int k = 0; k < i; k++) {
                    u[i * n + j] -= l[i * n + k] * u[k * n + j];
                }
            }
            else {
                l[i * n + j] = A(i,j);
                for (int k = 0; k < j; k++) {
                    l[i * n + j] -= l[i * n + k] * u[k * n + j];
                }
                l[i * n + j] = l[i * n + j]/u[j * n + j];
            }
        }

    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (lvl[i * n + j] == 0) {
                u[i * n + j] = 0;
                l[i * n + j] = 0;
            }
        }
    }


    Csr_matrix L(l);
    Csr_matrix U(u);
    L.Wolout(L);
    U.Wolout(U);
}

CSR_space::Csr_matrix CSR_space::Csr_matrix::CHOL0(const Csr_matrix &A) {
    int n = A.sup_row.size()-1;
    std::vector<double> l(n * n);
    std::vector<double> u(n * n);
    std::vector<double> lvl(n * n);


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (A(i,j) != 0) {
                lvl[i * n + j] = 1;
            }
        }
    }
//////////nice


    for (int i = 0; i <  n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += l[i * n + k] * l[j * n + k];

            if (i == j)
                l[i * n + j] = sqrt(A(i,i) - sum);
            else
                l[i * n + j] = (A(i,j) - sum) / l[j * n + j];
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u[i * n + j] = l[j * n + i];
        }
    }

    Csr_matrix L(l);
    Csr_matrix U(u);
    L.Wolout(L);
    U.Wolout(U);


    return L;
}

std::vector<double> CSR_space::Csr_matrix::SolverGaussReversed_forDownTriang(const Csr_matrix &A, const std::vector<double> &b) {
    int n = b.size();
    std::vector<double> x(n);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += A(i,j) * x[j];
        }
        x[i] = (b[i] - sum) / A(i,i);
    }
    return x;
}///working


std::vector<double> CSR_space::Csr_matrix::ReshConjugateGradientObusl(const Csr_matrix &A,const std::vector<double> &x01, const std::vector<double> &b,const double &e) {
    std::vector<double> x0 = x01;
    int n = A.sup_row.size() - 1;
    int iter = 0;

    double res = 1, alpha,beta;
    //need to do 1-st iteration
    std::vector<double> w0(n),r0(n), d0(n), rprev(n), wprev(n),dprev(n);
    r0 = A*x0 - b;
    Csr_matrix L = CHOL0(A);
    w0 = SolverGaussReversed_forDownTriang(L,r0);
    d0 = w0;

    while (res > e) {
        alpha = (r0 * w0) / Amult(A,d0);
        x0 = x0 - alpha * d0;
        rprev = r0;
        r0 = A* x0 - b;
        wprev = w0;
        w0 = SolverGaussReversed_forDownTriang(L,r0);
        res = vecmod(r0);
        beta = (r0 * w0) / (rprev * wprev);
        dprev = d0;
        d0 = w0 + beta * dprev;

        iter++;
    }
    std::cout << iter << "- number of iterations\n" ;
    return x0;
}




CSR_space::Csr_matrix::Csr_matrix(std::vector<double> vec1) {
    std::vector<DOK> vec;
    int n = std::sqrt(vec1.size());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vec.push_back({ i , j, vec1[i * n + j]});

        }
    }


    data.resize(vec.size());
    sup_row.resize(vec[vec.size() - 1].i + 2);
    col_ind.resize(vec.size());
    sup_row[0] = 0;

    std::size_t j = 1;
    for (std::size_t i = 0; i < vec.size(); i++) {
        data[i] = vec[i].val;
        col_ind[i] = vec[i].j;

        if (i > 0 && vec[i].i - vec[i - 1].i != 0) {
            for (int it = 0; it < vec[i].i - vec[i - 1].i; it++) {
                sup_row[j] = i;
                j++;
            }
        }
    }
    sup_row[j] = vec.size();
};
CSR_space::Csr_matrix::Csr_matrix(std::vector<DOK> vec) {
    data.resize(vec.size());
    sup_row.resize(vec[vec.size() - 1].i + 2);
    col_ind.resize(vec.size());
    sup_row[0] = 0;

    std::size_t j = 1;
    for (std::size_t i = 0; i < vec.size(); i++) {
        data[i] = vec[i].val;
        col_ind[i] = vec[i].j;

        if (i > 0 && vec[i].i - vec[i - 1].i != 0) {
            for (int it = 0; it < vec[i].i - vec[i - 1].i; it++) {
                sup_row[j] = i;
                j++;
            }
        }
    }
    sup_row[j] = vec.size();
};
