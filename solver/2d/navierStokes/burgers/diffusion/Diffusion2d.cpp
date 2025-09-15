#include <string>
#include <stdexcept>
#include "Diffusion2d.h"
#include "FieldUtil.h"

Diffusion2d::Diffusion2d(int meshX, int meshY, double constA, double deltaX, double deltaY, double deltaT) 
    : MESH_X(meshX), MESH_Y(meshY), CONST_A(constA), DELTA_X(deltaX), DELTA_Y(deltaY), DELTA_T(deltaT)
{
    validateTime();
}

void Diffusion2d::validateTime() {
    const auto DELTA = std::min(DELTA_X, DELTA_Y);
    const auto DELTA_T_MAX = 0.2 * DELTA * DELTA / CONST_A;
    if(DELTA_T > DELTA_T_MAX) {
        throw std::runtime_error("DELTA_T is too large. Must be <= " + std::to_string(DELTA_T_MAX));
    }
}

Field2d Diffusion2d::calculate(const Field2d& f) {
    Field2d f_next;
    FieldUtil::setSize(f_next, MESH_X, MESH_Y);

    for (int j = 1; j <= MESH_Y - 2; j++) {
        for (int i = 1; i <= MESH_X - 2; i++) {
            f_next[j][i] = f[j][i] + calculateTerm(f, i, j);
        }
    }
    return f_next;
}

Value Diffusion2d::calculateTerm(const Field2d& f, int i, int j) const {
    return (CONST_A * DELTA_T) *
           ( (f[j+1][i] - 2*f[j][i] + f[j-1][i]) / DELTA_Y / DELTA_Y + 
             (f[j][i+1] - 2*f[j][i] + f[j][i-1]) / DELTA_X / DELTA_X );
}
