#include <string>
#include <stdexcept>
#include "AdvectionDiffusion2d.h"
#include "FieldUtil.h"

AdvectionDiffusion2d::AdvectionDiffusion2d(int meshX, int meshY, double constA, double constU, double constV, double deltaX, double deltaY, double deltaT) 
    : MESH_X(meshX), MESH_Y(meshY), CONST_A(constA), CONST_U(constU), CONST_V(constV), DELTA_X(deltaX), DELTA_Y(deltaY), DELTA_T(deltaT),
      diffusion_(meshX, meshY, constA, deltaX, deltaY, deltaT), advection_(meshX, meshY, constU, constV, deltaX, deltaY, deltaT)
{
    validateTime();
}

void AdvectionDiffusion2d::validateTime() {
    diffusion_.validateTime();
    advection_.validateTime(CONST_U, CONST_V);
}

Field2d AdvectionDiffusion2d::calculate(const Field2d& f) {
    Field2d f_next;
    FieldUtil::setSize(f_next, MESH_X, MESH_Y);

    for (int j = 1; j <= MESH_Y - 2; j++) {
        for (int i = 1; i <= MESH_X - 2; i++) {
            f_next[j][i] = f[j][i] + advection_.calculateTerm(f, i, j) + diffusion_.calculateTerm(f, i, j);
        }
    }
    updateBoundaryCondition(f_next);
    return f_next;
}

void AdvectionDiffusion2d::updateBoundaryCondition(Field2d& f) {
    advection_.updateBoundaryCondition(f);
}