#include "Burgers2d.h"
#include <string>
#include <stdexcept>
#include "Advection2d.h"
#include "Diffusion2d.h"
#include "FieldUtil.h"

struct Burgers2d::Impl {
    Diffusion2d diffusion_;
    Advection2d advection_;

    Impl(int meshX, int meshY, double viscosity, double dx, double dy, double dt)
        : diffusion_(meshX, meshY, viscosity, dx, dy, dt),
          advection_(meshX, meshY, dx, dy, dt) {} 
};

Burgers2d::Burgers2d(int meshX, int meshY, double reynolds, 
                     double dx, double dy, double dt, Velocity2d f) 
    : MESH_X(meshX), MESH_Y(meshY), REYNOLDS(reynolds), 
      DX(dx), DY(dy), DT(dt), m_f(f),
      pimpl_(std::make_unique<Impl>(meshX, meshY, 1/reynolds, dx, dy, dt))
{
    validateTime();
}

Burgers2d::~Burgers2d() = default;

void Burgers2d::validateTime() {
    pimpl_->diffusion_.validateTime();

    const auto maxU = FieldUtil::findMax(m_f.u);
    const auto maxV = FieldUtil::findMax(m_f.v);
    pimpl_->advection_.validateTime(maxU, maxV);
}

Velocity2d Burgers2d::calculate() {
    Velocity2d f_next;
    FieldUtil::setSize(f_next.u, MESH_X, MESH_Y);
    FieldUtil::setSize(f_next.v, MESH_X, MESH_Y);

    for (int j = 1; j <= MESH_Y - 2; j++) {
        for (int i = 1; i <= MESH_X - 2; i++) {
            f_next.u[j][i] = m_f.u[j][i] + calculateTerm(m_f.u, m_f, i, j);
            f_next.v[j][i] = m_f.v[j][i] + calculateTerm(m_f.v, m_f, i, j);
        }
    }
    updateBoundaryCondition(f_next.u);
    updateBoundaryCondition(f_next.v);
    updateVelocity(f_next);
    return m_f;
}

Value Burgers2d::calculateTerm(const Field2d& f, const Velocity2d& velocity, int i, int j) const {
    return pimpl_->advection_.calculateVelocity(f, velocity.u, velocity.v, i, j) + pimpl_->diffusion_.calculateTerm(f, i, j);
}

void Burgers2d::updateBoundaryCondition(Field2d& f) {
    pimpl_->advection_.updateBoundaryCondition(f);
}

void Burgers2d::updateVelocity(Velocity2d& f) {
    std::swap(m_f, f);
}
