#pragma once
#include <memory>
#include "typedef.h"

class Burgers2d
{
public:
    Burgers2d(int meshX, int meshY, double reynolds, 
              double dx, double dy, double deltaT, Velocity2d f);
    ~Burgers2d();

    Velocity2d calculate();
    Value calculateTerm(const Field2d& f, const Velocity2d& velocity, int i, int j) const;
    void validateTime();
    void updateBoundaryCondition(Field2d& f);
private:
    void updateVelocity(Velocity2d& f);
    const int MESH_X, MESH_Y;
    const double REYNOLDS, DX, DY, DT;
    Velocity2d m_f;
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};
