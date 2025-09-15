#pragma once
#include "typedef.h"

class Poisson2d
{
public:
    Poisson2d(int meshX, int meshY, 
              double dx, double dy, double omega, double epsilon, double pRef, const MeshRange2d& range);
    ~Poisson2d() = default;

    int calculate(Field2d& p, const Field2d& s, const FlagField& ip, int iteration);
    Value calculateTerm(Field2d& p, const Field2d& s, const FlagField& ip) const;
private:
    const int MESH_X;
    const int MESH_Y;
    const double DX;
    const double DY;
    const double OMEGA;
    const double EPSILON;
    const double P_REF;
    const MeshRange2d MESH_RANGE;
};
