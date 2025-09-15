#include <string>
#include <stdexcept>
#include "Poisson2d.h"
#include "FieldUtil.h"

Poisson2d::Poisson2d(int meshX, int meshY, double dx, double dy, 
    double omega, double epsilon, double pRef, const MeshRange2d& range) 
    : MESH_X(meshX), MESH_Y(meshY), DX(dx), DY(dy), 
    OMEGA(omega), EPSILON(epsilon), P_REF(pRef), MESH_RANGE(range) 
{
    //do nothing
}


int Poisson2d::calculate(Field2d& p, const Field2d& s, const FlagField& ip, int iteration) {
    FieldUtil::ClearField(p);
    for (int time = 1; time <= iteration; time++) {
        const auto residual = calculateTerm(p, s, ip);
        if (residual < EPSILON ){
            return time;
        }
    }
    throw std::runtime_error("Poisson equation calculation results did not converge!!\nRe-consider poisson equation or Use higher iteration value!!");
    return iteration;
}

Value Poisson2d::calculateTerm(Field2d& p, const Field2d& s, const FlagField& ip) const {
    Value maxResidual = 0;
    for (int j = MESH_RANGE.minY; j <= MESH_RANGE.maxY; j++) {
        for (int i = MESH_RANGE.minX; i <= MESH_RANGE.maxX; i++) {
            if (ip[j][i] != ObjectFlag::fluid) {
                continue;
            }

            const auto p_n = (1 - OMEGA) * p[j][i] + OMEGA * 
            ((p[j][i + 1] + p[j][i - 1]) / DX / DX + (p[j + 1][i] + p[j - 1][i]) / DY / DY - s[j][i]) *
            DX * DX * DY * DY / (2 * (DX * DX + DY * DY));
            const auto currentResidual = std::abs(p_n - p[j][i]) / P_REF;
            maxResidual = std::max(maxResidual, currentResidual);
            p[j][i] = p_n;
        }
    }
    return maxResidual;
}
