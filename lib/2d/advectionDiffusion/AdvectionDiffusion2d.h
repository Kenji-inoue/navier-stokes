#pragma once
#include "typedef.h"
#include "Advection2d.h"
#include "Diffusion2d.h"


class AdvectionDiffusion2d
{
public:
    AdvectionDiffusion2d(int meshX, int meshY, double constA, double constU, double constV, double deltaX, double deltaY, double deltaT);
    ~AdvectionDiffusion2d() = default;

    Field2d calculate(const Field2d& f);
    void validateTime();
    void updateBoundaryCondition(Field2d& f);
private:
    const int MESH_X;
    const int MESH_Y;
    const double CONST_A;
    const double CONST_U;
    const double CONST_V;
    const double DELTA_X;
    const double DELTA_Y;
    const double DELTA_T;
    Diffusion2d diffusion_;
    Advection2d advection_;
};
