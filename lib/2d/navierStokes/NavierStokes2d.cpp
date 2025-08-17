#include <string>
#include <stdexcept>
#include "NavierStokes2d.h"
#include "FieldUtil.h"

NavierStokes2d::NavierStokes2d(int meshX, int meshY, double reynolds, 
                        double dx, double dy, double dt, double omega, double epsilon, double pRef, int poissonIteration,
                        const MeshRange2d& range, AnalysisResult& result, const Object& object)
    : MESH_X(meshX), MESH_Y(meshY), REYNOLDS(reynolds),
      DX(dx), DY(dy), DT(dt),
      EPSILON(epsilon), P_REF(pRef), OMEGA(omega), POISSON_ITERATION(poissonIteration), MESH_RANGE(range), m_result(result),
      burgers_(meshX, meshY, reynolds, dx, dy, DT, result.f),
      poisson_(meshX, meshY, dx, dy, omega, epsilon, pRef, range),
      m_fNext(result.f), m_object(object)
{
    m_dp.resize(MESH_Y, std::vector<Value>(MESH_X, 0.0));
}

AnalysisResult NavierStokes2d::calculate() {
    calculateProvisionalVelocity(m_fNext, m_result.f, m_result.p);
    calculateDivergenceOfVelocity(m_result.s, m_fNext);
    const auto interval = poisson_.calculate(m_dp, m_result.s, m_object.ip, POISSON_ITERATION);
    modifyPressure(m_result.p, m_dp);
    modifyVelocity(m_fNext, m_dp);
    calculateVorticity(m_result.rot, m_fNext);
    calculateDragForce(m_result.drag, m_fNext, m_result.p);

    updateVelocityTimeScale();
    return m_result;
}

void NavierStokes2d::calculateProvisionalVelocity(Velocity2d& fNext, const Velocity2d& f, const Field2d& p) {
    for (int j = MESH_RANGE.minY+1; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX+1; i <= MESH_RANGE.maxX; i++) {
            if (m_object.iu[j][i] != ObjectFlag::fluid) {
                continue;
            }

            fNext.u[j][i] = f.u[j][i] + burgers_.calculateTerm(f.u, f, i, j) - calculatePressureTermX(p, i, j);
            fNext.v[j][i] = f.v[j][i] + burgers_.calculateTerm(f.v, f, i, j) - calculatePressureTermY(p, i, j);
        }
    }
    updateRunoffBoundaryCondition(fNext);
}

Value NavierStokes2d::calculatePressureTermX(const Field2d& p, int i, int j) const {
    const Value frontPressure = (p[j][i] + p[j - 1][i]) / 2;
    const Value backPressure = (p[j][i - 1] + p[j - 1][i - 1]) / 2;
    return DT * (frontPressure - backPressure) / DX;
}

Value NavierStokes2d::calculatePressureTermY(const Field2d& p, int i, int j) const {
    const Value frontPressure = (p[j][i] + p[j][i - 1]) / 2;
    const Value backPressure = (p[j - 1][i] + p[j - 1][i - 1]) / 2;
    return DT * (frontPressure - backPressure) / DY;
}

void NavierStokes2d::calculateDivergenceOfVelocity(Field2d& s, const Velocity2d& f) {
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX; i <= MESH_RANGE.maxX; i++) {
            if (m_object.ip[j][i] != ObjectFlag::fluid) {
                continue;
            }

            s[j][i] = (((f.u[j + 1][i + 1] - f.u[j + 1][i]) / DX + (f.u[j][i + 1] - f.u[j][i]) / DX) / 2 +
                      ((f.v[j + 1][i + 1] - f.v[j][i + 1]) / DY + (f.v[j + 1][i] - f.v[j][i]) / DY) / 2) / DT;
        }
    }
}

void NavierStokes2d::updateRunoffBoundaryCondition(Velocity2d& f) {
    for (int j = MESH_RANGE.minY+1; j <= MESH_RANGE.maxY; j++) {
        f.u[j][MESH_X - 2] = f.u[j][MESH_X - 3];
        f.u[j][MESH_X - 1] = f.u[j][MESH_X - 3];
        f.v[j][MESH_X - 2] = f.v[j][MESH_X - 3];
        f.v[j][MESH_X - 1] = f.v[j][MESH_X - 3];
    }
}

void NavierStokes2d::modifyPressure(Field2d& p, Field2d& dp) {
    for (int j = 0; j < MESH_Y; j++) {
        for (int i = 0; i < MESH_X; i++) {
            if (m_object.ip[j][i] != ObjectFlag::fluid) {
                continue;
            }

            p[j][i] += dp[j][i];
        }
    }

    // Runoff Boundary Condition
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        dp[j][MESH_X - 2] = p[j][MESH_X - 3] - p[j][MESH_X - 2];
        p[j][MESH_X - 2] = p[j][MESH_X - 3];
    }

    // Inflow Boundary Condition
    for (int j = 0; j < MESH_Y; j++) {
        p[j][0] = 0;
        p[j][1] = 0;
    }
    for (int i = 0; i < MESH_X; i++) {
        p[0][i] = 0;
        p[1][i] = 0;
        p[MESH_Y - 1][i] = 0;
        p[MESH_Y - 2][i] = 0;
    }
}   

void NavierStokes2d::modifyVelocity(Velocity2d& f, const Field2d& dp) {
    for (int j = MESH_RANGE.minY+1; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX+1; i <= MESH_RANGE.maxX; i++) {
            if (m_object.iu[j][i] != ObjectFlag::fluid) {
                continue;
            }

            f.u[j][i] = f.u[j][i] - DT / 2 * ((dp[j][i] - dp[j][i - 1]) / DX + (dp[j - 1][i] - dp[j - 1][i - 1]) / DX);
            f.v[j][i] = f.v[j][i] - DT / 2 * ((dp[j][i] - dp[j - 1][i]) / DY + (dp[j][i - 1] - dp[j - 1][i - 1]) / DY);
        }
    }
    updateRunoffBoundaryCondition(f);
}

void NavierStokes2d::calculateVorticity(Field2d& rot, const Velocity2d& f) {
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX; i <= MESH_RANGE.maxX; i++) {
            if (m_object.ip[j][i] != ObjectFlag::fluid) {
                continue;
            }

            rot[j][i] = ((f.v[j + 1][i + 1] - f.v[j + 1][i] + f.v[j][i + 1] - f.v[j][i]) / DX -
                         (f.u[j + 1][i + 1] - f.u[j][i + 1] + f.u[j + 1][i] - f.u[j][i]) / DY) / 2;
        }
    }
}

void NavierStokes2d::calculateDragForce(DragForce& drag, const Velocity2d& f, const Field2d& p) {
    drag.x = 0.0;
    drag.y = 0.0;
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX; i <= MESH_RANGE.maxX; i++) {
            if (m_object.ip[j][i] != ObjectFlag::fluid) {
                continue;
            }

            const auto dudx = (f.u[j + 1][i + 1] - f.u[j + 1][i] + f.u[j][i + 1] - f.u[j][i]) / 2.0 / DX;
            const auto dudy = (f.u[j + 1][i + 1] + f.u[j + 1][i] - f.u[j][i + 1] - f.u[j][i]) / 2.0 / DY;
            const auto dvdx = (f.v[j + 1][i + 1] - f.v[j + 1][i] + f.v[j][i + 1] - f.v[j][i]) / 2.0 / DX;
            const auto dvdy = (f.v[j + 1][i + 1] + f.v[j + 1][i] - f.v[j][i + 1] - f.v[j][i]) / 2.0 / DY;

            const auto sigmaX = -p[j][i] + 4.0 / 3.0 * dudx / REYNOLDS - 2.0 / 3.0 * dvdy / REYNOLDS;
            const auto tauX = (dvdx + dudy) / REYNOLDS;
            const auto sigmaY = -p[j][i] + 4.0 / 3.0 * dvdy / REYNOLDS - 2.0 / 3.0 * dudx / REYNOLDS;
            const auto tauY = tauX;

            if (m_object.ip[j][i + 1] == ObjectFlag::surface) {
                drag.x -= sigmaX * DY;
                drag.y -= tauX * DY;
            }
            else if (m_object.ip[j][i - 1] == ObjectFlag::surface) {
                drag.x += sigmaX * DY;
                drag.y += tauX * DY;
            }

            if (m_object.ip[j + 1][i] == ObjectFlag::surface) {
                drag.x -= tauY * DX;
                drag.y -= sigmaY * DX;
            }
            else if (m_object.ip[j - 1][i] == ObjectFlag::surface) {
                drag.x += tauY * DX;
                drag.y += sigmaY * DX;
            }
        }
    }
}

void NavierStokes2d::updateVelocityTimeScale() {
    std::swap(m_result.f, m_fNext);
}
