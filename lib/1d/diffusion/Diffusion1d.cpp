#include <string>
#include <stdexcept>
#include "Diffusion1d.h"

Diffusion1d::Diffusion1d(int meshX, double constA, double deltaX, double deltaT) 
    : MESH_X(meshX), CONST_A(constA), DELTA_X(deltaX), DELTA_T(deltaT)
{
    validateTime();
}

void Diffusion1d::validateTime() {
    const auto DELTA_T_MAX = 0.2 * DELTA_X * DELTA_X / CONST_A;
    if(DELTA_T > DELTA_T_MAX) {
        throw std::runtime_error("DELTA_T is too large. Must be <= " + std::to_string(DELTA_T_MAX));
    }
}

Field1d Diffusion1d::calculate(const Field1d& f) {
    Field1d f_next(MESH_X);
    const auto coeff = (CONST_A * DELTA_T / DELTA_X / DELTA_X);
    for (int i = 1; i <= MESH_X - 2; i++) {
        f_next[i] = f[i] + coeff * (f[i + 1] - 2 * f[i] + f[i - 1]);
    }
    return f_next;
}
