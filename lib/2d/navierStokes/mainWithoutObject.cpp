#include "NavierStokes2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void setInflowBoundaryCondition(Field2d& f, int meshX, int meshY) {
    for (int j = 0; j < meshY; j++) {
        f[j][0] = 1.0;
        f[j][1] = 1.0;
    }
}

void defineObject(Object& object, int meshX, int meshY) {
    FieldUtil::InitializeFlagField(object.iu, meshX, meshY, ObjectFlag::fluid);
    FieldUtil::InitializeFlagField(object.ip, meshX, meshY, ObjectFlag::fluid);
}

int main() {
    int meshX = 16 + 3;
    int meshY = 8 + 3;
    double reynolds = 200;
    double dx = 1.0 / (meshX - 3);
    double dy = dx;
    double dt = 0.2 * dx / 1.0;
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
    setInflowBoundaryCondition(result.f.u, meshX, meshY);

    FieldUtil::InitializeField(result.rot, meshX, meshY, 0);
    result.drag.x = 0.0;
    result.drag.y = 0.0;

    Object object;
    defineObject(object, meshX, meshY);

    try {
        NavierStokes2d solver(meshX, meshY, reynolds, dx, dy, dt,
                            omega, epsilon, pRef, poissonIteration, 
                            range, result, object);
        const int interval = 1;
        const int maxIterations = 2;
        for (int time = 0; time < maxIterations; time++) {
            printf("Iteration: %d\n", time);
            printf("u result\n");
            FieldUtil::display(result.f.u, time, interval);

            printf("v result\n");
            FieldUtil::display(result.f.v, time, interval);

            printf("p result\n");
            FieldUtil::display(result.p, time, interval);

            printf("s result\n");
            FieldUtil::display(result.s, time, interval);
            result = solver.calculate();
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}