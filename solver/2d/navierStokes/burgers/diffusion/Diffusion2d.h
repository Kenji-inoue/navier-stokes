#pragma once
#include "typedef.h"

class Diffusion2d
{
public:
    Diffusion2d(int meshX, int meshY, double constA, double deltaX, double deltaY, double deltaT);
    ~Diffusion2d() = default;

    Field2d calculate(const Field2d& f);

    Value calculateTerm(const Field2d& f, int i, int j) const;
    void validateTime();
private:
    const int MESH_X;
    const int MESH_Y;
    const double CONST_A;
    const double DELTA_X;
    const double DELTA_Y;
    const double DELTA_T;
};
