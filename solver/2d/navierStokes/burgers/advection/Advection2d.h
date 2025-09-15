#pragma once
#include "typedef.h"

class Advection2d
{
public:
    Advection2d(int meshX, int meshY, double constU, double constV, double deltaX, double deltaY, double deltaT);
    Advection2d(int meshX, int meshY, double deltaX, double deltaY, double deltaT);
    ~Advection2d() = default;

    Field2d calculate(const Field2d& f);

    Value calculateTerm(const Field2d& f, int i, int j) const;
    Value calculateVelocity(const Field2d& f, const Field2d& u, const Field2d& v, int i, int j) const;
    void validateTime(Value u, Value v);
    void updateBoundaryCondition(Field2d& f);
private:
    const int MESH_X;
    const int MESH_Y;
    const double CONST_U;
    const double CONST_V;
    const double DELTA_X;
    const double DELTA_Y;
    const double DELTA_T;
};
