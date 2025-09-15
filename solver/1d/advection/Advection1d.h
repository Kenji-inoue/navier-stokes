#pragma once
#include "typedef.h"

class Advection1d
{
public:
    Advection1d(int meshX, double constA, double deltaX, double deltaT);
    ~Advection1d() = default;

    Field1d calculate(const Field1d& f);
private:
    void validateTime();
    void periodicBoundaryCondition(Field1d& f);
    const int MESH_X;
    const double CONST_U;
    const double DELTA_X;
    const double DELTA_T;
};
