//
// Created by gorilla on 10.02.23.
//

#include "tridiagonal.h"

void TriMat::apply(const std::vector<line>& data1) {
 data = data1;

}

const line& TriMat::operator()(int i) const{
    return data[i];
}


line& TriMat::operator()(int i) {
    return data[i];
}
