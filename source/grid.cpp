#include <vector>
#include "../include/grid.hpp"

size_t grid3D::getIndex(int i, int j, int k) {
    return size_t(k + this->N.z*(j + this->N.y*i));
}

grid3D::grid3D(pod3<int> N, pod3<double> L) : F(N.x*N.y*N.z) {
    this->N = N;
    this->L = L;
}

void grid3D::set(int i, int j, int k, double val) {
    size_t index = grid3D::getIndex(i, j, k);
    this->F[index] = val;
}

double grid3D::get(int i, int j, int k) {
    size_t index = grid3D::getIndex(i, j, k);
    return this->F[index];
}
