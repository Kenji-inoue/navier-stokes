#pragma once
#include "typedef.h"
#include "Burgers2d.h"
#include "Poisson2d.h"

struct DragForce {
    double x;
    double y;
};

struct AnalysisResult {
    Velocity2d f;
    Field2d p;
    Field2d s;
    Field2d rot;
    DragForce drag;
};

class NavierStokes2d
{
public:
    NavierStokes2d(int meshX, int meshY, double reynolds, double dx, double dy, double dt,
                double omega, double epsilon, double pRef, int poissonIteration,
                const MeshRange2d& range, AnalysisResult& result, const Object& object);
    ~NavierStokes2d() = default;
    AnalysisResult calculate();

private:
    void calculateProvisionalVelocity(Velocity2d& fNext, const Velocity2d& f, const Field2d& p);
    Value calculatePressureTermX(const Field2d& p, int i, int j) const;
    Value calculatePressureTermY(const Field2d& p, int i, int j) const;
    void calculateDivergenceOfVelocity(Field2d& s, const Velocity2d& f);
    void updateRunoffBoundaryCondition(Velocity2d& f);
    void modifyPressure(Field2d& p, Field2d& dp);
    void modifyVelocity(Velocity2d& f, const Field2d& dp);
    void calculateVorticity(Field2d& rot, const Velocity2d& f);
    void calculateDragForce(DragForce& drag, const Velocity2d& f, const Field2d& p);
    void updateVelocityTimeScale();
    const int MESH_X;
    const int MESH_Y;
    const double REYNOLDS;
    const double DX;
    const double DY;
    const double DT;
    const double EPSILON;
    const double P_REF;
    const double OMEGA;
    const int POISSON_ITERATION;
    const MeshRange2d MESH_RANGE;
    AnalysisResult m_result;
    Velocity2d m_fNext;
    Field2d m_dp;
    Burgers2d burgers_;
    Poisson2d poisson_;
    Object m_object;
};
