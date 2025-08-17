#include "NavierStokes2d.h"
#include "FieldUtil.h"
#include "FileUtil.h"
#include <stdexcept>
#include <iostream>

void setInflowBoundaryCondition(Velocity2d& f, int meshX, int meshY) {
    for (int j = 0; j < meshY; j++) {
        f.u[j][0] = 0.98;
        f.u[j][1] = 0.98;
        f.v[j][0] = 0.02;
        f.v[j][1] = 0.02;
    }
    for (int i = 0; i < meshX; i++) {
        f.u[0][i] = 0.98;
        f.u[1][i] = 0.98;
        f.u[meshY - 1][i] = 0.98;
        f.u[meshY - 2][i] = 0.98;

        f.v[0][i] = 0.02;
        f.v[1][i] = 0.02;
        f.v[meshY - 1][i] = 0.02;
        f.v[meshY - 2][i] = 0.02;
    }
}

void defineObject(Object& object, int meshX, int meshY) {
    FieldUtil::InitializeFlagField(object.iu, meshX, meshY, ObjectFlag::fluid);
    FieldUtil::InitializeFlagField(object.ip, meshX, meshY, ObjectFlag::fluid);
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            if (28 <= j && j <= 36 && 28 <= i && i <= 36) {
                object.iu[j][i] = ObjectFlag::surface;
            }
            if (28 < j && j < 36 && 28 < i && i < 36) {
                object.iu[j][i] = ObjectFlag::inside;
            }
            if (28 <= j && j <= 35 && 28 <= i && i <= 35) {
                object.ip[j][i] = ObjectFlag::surface;
            }
            if (28 < j && j < 35 && 28 < i && i < 35) {
                object.ip[j][i] = ObjectFlag::inside;
            }
        }
    }
}

int main() {
    int meshX = 128 + 3;
    int meshY = 64 + 3;
    double reynolds = 200;
    double dx = 16.0 / (Value)(meshX - 3);
    double dy = dx;
    double dt = 0.5 * dx / 4.0;
    AnalysisResult result;
    double omega = 1.0;
    double epsilon = 1e-7;
    double pRef = 1.0;
    int poissonIteration = 99999;
    MeshRange2d range = {1, meshX - 3, 1, meshY - 3};
    FieldUtil::InitializeField(result.f.u, meshX, meshY, 0);
    FieldUtil::InitializeField(result.f.v, meshX, meshY, 0);
    FieldUtil::InitializeField(result.p, meshX, meshY, 0);
    FieldUtil::InitializeField(result.s, meshX, meshY, 0);
    FieldUtil::InitializeField(result.rot, meshX, meshY, 0);
    setInflowBoundaryCondition(result.f, meshX, meshY);
    result.drag.x = 0.0;
    result.drag.y = 0.0;

    Object object;
    defineObject(object, meshX, meshY);

    try {
        NavierStokes2d solver(meshX, meshY, reynolds, dx, dy, dt,
                             omega, epsilon, pRef, poissonIteration, 
                             range, result, object);
        const int interval = 50;
        const int maxIterations = 500;
        FileUtil file("result.csv");
        for (int time = 1; time <= maxIterations; time++) {
            result = solver.calculate();

            if (time % interval != 0) {
                continue;
            }

            file.saveField(result.rot, "rot", time, interval);
            printf("dragX = %6.3f, dragY = %6.3f\n", result.drag.x, result.drag.y);
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}