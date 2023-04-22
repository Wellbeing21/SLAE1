//
// Created by gorilla on 10.02.23.
//

#ifndef SLAE_TRIDIAGONAL_H
#define SLAE_TRIDIAGONAL_H
#include <vector>
#include <array>
struct line {
    double f;
    double s;
    double t;
};
class TriMat {
private:
    //private field
    std::vector <line> data;
public: //pity
    void apply(const std::vector<line>& data1);
    //overdosing operator
    const line& operator()(int i) const;
    line& operator() (int i);


};




#endif //SLAE_TRIDIAGONAL_H
