#include "dense.h"

/////////////////////
double vecmod (const std::vector<double> &a) {
    double s = 0;
    for (int i = 0; i < a.size(); i++) {
        s += a[i] * a[i];
    }
    return sqrt(s);
}

void out(std::vector<double> &a) {
    int n = ( a.size());
    std::cout << "vec:\n";
    for (int i = 0; i < n; i++) {
        std::cout << a[i]<<"\n";
    }
    std::cout << "\n";
}

double operator* (const std::vector<double> &a, const std::vector<double> &b) {
    double sum = 0;
    for (int i = 0; i < a.size(); i++) {
        sum += a[i] * b[i];
    }
    return sum;
}
/////////////////////


    Matr::Matr(std::vector<std::vector<double>> &a) {
        //with vector of vectors
        int n = a.size();
        vv.resize(n * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vv[i * n + j] = a[i][j];
            }
        }
    }
    Matr::Matr(int n) {
        //with zeros
        vv.resize(n * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vv[i * n + j] = 0;
            }
        }
    }
    Matr::Matr (std::vector <double> &a) {
        int n = a.size();
        vv.resize(n);
        for (int i = 0; i < n; i++) {
            vv[i] = a[i];
        }
    }
    void Matr::out() {
        int n = sqrt( vv.size());
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << std::setw(10) << vv[i * n + j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    double Matr::operator()(const int &i,const int &j) const{
        int n = std::sqrt(vv.size());
        return vv[i * n + j];
    }
//with const if we don't need to change something
    double Matr::operator()(const int &i,const  int &j) {
        int n = sqrt(vv.size());
        return vv[i * n + j];
    }
//without(to change smth)
    Matr Matr::operator* (Matr &M) {
        int n = sqrt(M.vv.size());
        Matr m(n);
//created empty matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    m.vv[i * n + j] += vv[i * n + k] * M.vv[k * n + j];
                }
            }
        }
        return m;
    }
    std::pair<Matr, Matr> Matr::QRH () {
        int n = sqrt(vv.size());
        Matr Q(n);
        //try to make Q matrix
        double modm = 0;
        std::vector<double> x(n);
        for (int i = 0; i < n; i++){
            x[i] = vv[i * n];
            modm += vv[i * n] * vv[i * n];
        }
        modm = sqrt(modm);
        x[0] += modm;
        modm = 0;
        //made 1-st vector normal
        for (int i = 0; i < n; i++) {
            Q.vv[i * n + i] += 1;
            modm += x[i] * x[i];
        }
        //made I and |V|^2
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q.vv[i * n + j] -= 2 / modm * x[i] * x[j];
            }
        }
        //YEAH i made P1

        for (int ht = 0; ht < n - 1; ht++) {
            //householder times

            std::vector<double> nor(n);
            double mod = 0;
            for (int i = ht; i < n; i++) {
                mod += vv[i * n + ht] * vv[i * n + ht];
            }
            mod = sqrt(mod);
            //take our |x|
            for (int i = ht; i < n; i++) {
                nor[i] = vv[i * n + ht];
            }
            nor[ht] += mod;
            //created normal vector
            mod = 0;
            for (int i = ht; i < n; i++) {
                mod += nor[i] * nor[i];
            }
            //created |V|^2
            ////was for (int vn = ht; vn < n; vn++)
            for (int vn = 0; vn < n; vn++) {
                //vector number
                double sc = 0;
                //scalar multiplication
                for (int i = ht; i < n; i++) {
                    sc += vv[i * n + vn] * nor[i];
                }
                //finded our scalar multiplication for every colomn vector
                for (int i = ht; i < n; i++) {
                    vv[i * n + vn] -= 2 * sc / mod * nor[i];
                    // was Q.vv[i  * n + vn] -= 2 * sc / mod * nor[i];
                }
                sc = 0;
                if (ht > 0) {
                    for (int i = ht; i < n; i++) {
                        sc += Q.vv[vn * n + i] * nor[i];
                    }

                    for (int i = ht; i < n; i++) {
                        Q.vv[vn  * n + i] -= 2 * sc / mod * nor[i];
                        //but i can't do it like this, because scalar multiplication isn't the same
                    }
                }

            }
        }
        //and we need to transpose
        Matr R(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ( i > j && fabs(R.vv[i * n + j]) < 0.00000000001) {
                    R.vv[i * n + j] = 0;
                }
                else {R.vv[i * n + j] =  vv[i * n + j];}
            }
        }
        return std::make_pair(Q, R);


    }
    std::pair<Matr, Matr> Matr::QRG () {
        int n = sqrt(vv.size());
        //here we need to check if [0,0] ne 0, so
        if (vv[0] == 0) {
            bool stop = false;
            int i = 0;
            while (vv[i * n] == 0) {
                if (stop) {
                    break;
                }
                i++;
                if (vv[i * n] != 0) {
                    stop = true;
                    //take givence rotation with this element
                    for (int j = 0; j < n; j++) {
                        double temp;
                        temp = vv[0 * n + j];
                        vv[0 * n + j] = vv[i * n + j];
                        vv[i * n + j] = -temp;

                    }
                }
            }

        }
        ////its working: if we have zeros in first colomn, they will rotate to value
        Matr Q(n);
        for (int i = 0; i < n; i++) {
            Q.vv[i * n + i] = 1;
        }
        //created P1
        double c = vv[0] / sqrt(vv[0] * vv[0] + vv[1 * n] * vv[1 * n]), s = vv[1 * n] / sqrt(vv[0] * vv[0] + vv[1 * n] * vv[1 * n]);
        Q.vv[0] = c;
        Q.vv[1 * n + 1] = c;
        Q.vv[0 * n + 1] = s;
        Q.vv[1 * n + 0] = -s;
        ////P1's working
        for (int gt = 0; gt < n; gt++) {
            for (int sn = gt + 1; sn < n; sn++) {

                //here we have Givence time and string number
                double sinn = vv[sn * n + gt] / sqrt(vv[sn * n + gt] * vv[sn * n + gt] + vv[gt * n + gt] * vv[gt * n + gt]);
                double coss = vv[gt * n + gt] / sqrt(vv[sn * n + gt] * vv[sn * n + gt] + vv[gt * n + gt] * vv[gt * n + gt]);
                //for every rotation we need to find our sinn and coss(for every string)
                for (int j = 0; j < n; j++) {
                    double base = 0, string = 0, baseq = 0, stringq = 0;
                    base += vv[gt * n + j] * coss + vv[sn * n + j] * sinn;
                    string += vv[sn * n + j] * coss - vv[gt * n + j] * sinn;
                    //temporarily to multiply
                    //i - nomer coordinatu
                    //starts with next
                    vv[gt * n + j] = base;
                    vv[sn * n + j] = string;
                    if ((sn > gt + 1) or (gt > 0)) {
                        baseq += Q.vv[gt * n + j] * coss + Q.vv[sn * n + j] * sinn;
                        stringq += Q.vv[sn * n + j] * coss - Q.vv[gt * n + j] * sinn;
                        Q.vv[gt * n + j] = baseq;
                        Q.vv[sn * n + j] = stringq;
                    }

                    //transform our sn string here
                }

                //we multiplied here our sn string on matrix, but also first string
////its working also


            }
        }
        Matr R(vv);
        return std::make_pair(Q, R);
    }




std::vector<double>  Matr::operator* (const std::vector<double> &a) const {
    int n = a.size();
    std::vector<double> x0(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            x0[i] += vv[i * n + j] * a[j];
        }
    }
    return x0;

}




std::vector<double> Matr::SolverGaussReversed_forUpTriang(const Matr &A, const std::vector<double> &b) {
    int n = b.size();
    std::vector<double> x(n);
    for (int i = n-1; i >= 0; i--) {
        double sum = 0;
        for (int j = i+1; j < n; j++) {
            sum += A(i,j) * x[j];
        }
        x[i] = (b[i] - sum) / A(i,i);
    }
    return x;
}///working
std::vector<double> Matr::GMRES ( const Matr &A,  const std::vector<double> &x01, const std::vector<double> &b,const double &e,const int &m1) {
    int n = b.size();
    int it = 0;
    int jt = 0;
    double ri = 10;
    std::vector<double> zi(m1);
    std::vector<double> x0 = x01;
    double m = m1 - 1;
    while (not (ri < e)) {
        double t = 0, it = 0;
        std::vector<double> Givens(2 * m);
        std::vector<std::vector<double>> vi(n);//? maybe not n * n
        std::vector<double> H((m+1) * m);//say my name, Heisenberg, you god damn right
        for (int i = 0; i < n; i++) {
            vi[i].resize(n);
        }
        std::vector<double> r = A * x0;
        for (int i = 0; i < n; i++) {
            r[i] = r[i] - b[i];
        }
        double modr = vecmod(r);
        double modr0 = modr;
        zi[0] = modr0;
        for (int i = 0; i < n; i++) {
            r[i] = r[i] / modr;
        }
        vi[0] = r;
        //created v0
        for (int k = 1; k < m + 1; k++) {
            vi[k] = A * vi[k - 1];
            //next
            for (int i = 0; i < k; i++) {
                H[i * m + (k - 1)] = (vi[k]*vi[i]);
            }
            //filled all include i, but not delta i+1, i
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < n; j++) {
                    vi[k][j] -= H[i * m + (k - 1)] * vi[i][j];
                }
            }
            modr = vecmod(vi[k]);
            H[k * m + (k - 1)] = modr;
            for (int i = 0; i < n; i++) {
                vi[k][i] = vi[k][i] / modr;
            }
            //tut nado yzhe sdelat vrasheniya
            it = 0;
            for (int i = 0; i < k-1; i++) {
                double a = H[(i) * m + k-1] * Givens[it] + H[(i+1) * m + k-1] * Givens[it + 1];
                double b = -H[(i) * m + k-1] * Givens[it + 1] + H[(i+1) * m + k-1] * Givens[it];
                H[(i) * m + k-1] = a;
                H[(i+1) * m + k-1] = b;
                it += 2;
            }


            double squareroot = sqrt(
                    H[t * m + k - 1] * H[t * m + k - 1] + H[(t + 1) * m + k - 1] * H[(t + 1) * m + k - 1]);
            Givens[it] = H[t * m + k - 1] / squareroot; //cos delta00/squareroot
            t++;
            it++;
            Givens[it] = H[t * m + k - 1] / squareroot; //sin delta 01/squareroot
            it++;
            //now we need to do that on Hessinberg matrix
            //final vrasheniye
            it-=2;
            double a = H[(k-1) * m + k-1] * Givens[it] + H[(k) * m + k-1] * Givens[it + 1];
            double b = -H[(k-1) * m + k-1] * Givens[it + 1] + H[(k) * m + k-1] * Givens[it];
            H[(k-1) * m + k-1] = a;
            H[(k) * m + k-1] = b;
        }
        //zapisali vse vrasheniya



        //sdelali eto
        //filled hessinberg matrix, but i DONT THOUGHT OF LAST Symbol

        //vse gotovo, we need to do givens and then find solution

        it =0;

        for (int i = 0; i < m; i++) {
            double a = zi[i] * Givens[it] + zi[i + 1] * Givens[it + 1];
            double b = -zi[i] * Givens[it + 1] + zi[i + 1] * Givens[it];
            it += 2;
            zi[i] = a;
            zi[i + 1] = b;
        }
        ri = fabs(zi[m]);
        //created zi

        //   H [m*m - 1] = H[m * m - 1] * Givens[2 * m - 2] + delta * Givens [2 * m - 1];
        Matr H1 = Matr(H);
        std::vector<double> zibezlast(m);
        for (int i = 0; i < m; i++) {
            zibezlast[i] = zi[i];
        }
        std::vector<double> y = H1.SolverGaussReversed_forUpTriang(H1,zibezlast);
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < n; k++) {
                x0[k] -= vi[i][k] * y[i];
            }
        }
        jt++;
    }

    return (x0);
}

std::vector<double> Matr::TransposedMatrixMult (const std::vector<double> &a) const {
    std::vector<double>result(a.size());
    for (int k = 0; k < a.size(); k++) {
        for(int i = 0; i < a.size(); i++) {
            result[k] += vv[i * a.size() + k] * a[i];
        }
    }



    return result;
}

std::vector<double> Matr::BCG ( const Matr &A,  const std::vector<double> &x0, const std::vector<double> &b,const double e) {
    double q, th, res = 10, n = x0.size();
    std::vector<double> ri2, xi = x0;
    ri2 = A*x0;
    for (int i = 0; i < n; i++) {
        ri2[i] -=  b[i];
    }
    double it = 0;
    //found r0
    std::vector<double> rti2 = ri2, rti1 = ri2, ri1 = ri2, pi = ri2, pti = ri2;
    while (fabs(res) > e) {
        it++;
        std::vector<double> Api = A*pi, Atpti = A.TransposedMatrixMult(pti);
        q = (rti1*ri1) / (rti1*( A*pi));
        for (int i = 0; i < n; i++) {
            xi[i] = xi[i] - q * pi[i];
        }
        //i'm going to recalculate residual so i write it to ri2
        ri2 = ri1;
        rti2 = rti1;
        for (int i = 0; i < n; i++) {
            ri1[i] = ri1[i] - q * Api[i];
            rti1[i] = rti1[i] - q * Atpti[i];
        }
        //found xi,ri,rti, now theta
        th = (rti1*ri1) / (rti2*ri2);
        // now computing pi, pti
        for (int i = 0; i < n; i++) {
            pi[i] = ri1[i] + th * pi[i];
            pti[i] = rti1[i] + th * pti[i];
        }
        res = vecmod(ri1);



    }
    std::cout << it << " iterations \n";
    return xi;
}


std::vector<double> Matr::CGS ( const Matr &A,  const std::vector<double> &x0, const std::vector<double> &b,const double e) {
    double q, th, res = 10, n = x0.size();
    std::vector<double> ris, xi = x0;
    double it = 0;
    while (fabs(res) > e) {
        it = 0;
        ris = A*xi;
        for (int i = 0; i < n; i++) {
            ris[i] -=  b[i];
        }
        //found r0
        std::vector<double> h(n), w(n), pis = ris, rt0 = ris, ui = ris, rsiprev(n);
        while (fabs(res) > e) {
            it++;
            if (it > 1) {
                th = (ris*ris) / (rsiprev*rsiprev);
                for (int i = 0; i < n; i++) {
                    ui[i] = ris[i] + th * h[i];
                    w[i] = q + th * pis[i];
                    pis[i] = ui[i] + th * w[i];
                }
            }
            w = A * pis;
            if ((ris*w) == 0) { q = 0; break;}
            q = (ris*ris) / (pis*w);
            for (int i = 0; i < n; i++) {
                h[i] = ui[i] - q * w[i];
                w[i] = ui[i] + h[i];
            }
            rsiprev = ris;
            std::vector<double> Aw = A * w;
            for (int i = 0; i < n; i++) {
                ris[i] = ris[i] - q * Aw[i];
                xi[i] = xi[i] - q * w[i];
            }
            res = vecmod(ris);
            if (fabs(q) < 0.000005) { break;}



        }
    }
    return xi;
}


