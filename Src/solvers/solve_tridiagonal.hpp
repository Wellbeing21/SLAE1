//
// Created by gorilla on 10.02.23.
//
#include "../Tridiagonal_matrix/tridiagonal.h"

namespace solver {
    std::vector<double> findsolution(const TriMat &Matr, std::vector<double> free) {
        int n = free.size();
        std::vector<double> p_(n);
        std::vector<double> q_(n);
        std::vector<double> sol(n);
        p_[0] = -Matr(0).t / Matr(0).s;
        q_[0] = -free[0] / Matr(0).s;
        for (int i = 1; i < n; i++) {
            p_[i] = -Matr(i - 1).t / (Matr(i - 1).f * p_[i - 1] + Matr(i - 1).s);
            q_[i] = (free[i - 1] - Matr(i - 1).f * q_[i - 1]) / (Matr(i - 1).f * p_[i - 1] + Matr(i - 1).s);
        }
        sol[n - 1] = (free[n - 1] - Matr(n - 1).f * q_[n - 1]) / (Matr(n - 1).f * p_[n - 1] + Matr(n - 1).s);
        for (int i = n - 2; i >= 0; i--) {
            sol[i] = p_[i + 1] * sol[i + 1] + q_[i + 1];

        }
        return sol;
    }
}


