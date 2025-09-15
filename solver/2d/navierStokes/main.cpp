#include "NavierStokes2d.h"
#include "FieldUtil.h"
#include "FileUtil.h"
#include "Hdf5Util.h"
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include <string>

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

void readConfig(Hdf5Util& file, int& meshX, int& meshY, Value& dx, Value& dy, Value& dt,
                    Value& reynolds, Value& omega, Value& epsilon, Value& pRef,
                    int& poissonIteration, int& maxIterations, int& interval) {
    file.readIntConfig("meshX", meshX);
    file.readIntConfig("meshY", meshY);
    file.readDoubleConfig("reynolds", reynolds);
    file.readDoubleConfig("dx", dx);
    file.readDoubleConfig("dy", dy);
    file.readDoubleConfig("dt", dt);
    file.readDoubleConfig("omega", omega);
    file.readDoubleConfig("epsilon", epsilon);
    file.readDoubleConfig("pRef", pRef);
    file.readIntConfig("poissonIteration", poissonIteration);
    file.readIntConfig("maxIterations", maxIterations);
    file.readIntConfig("interval", interval);
}

void InitializeAnalysisResult(AnalysisResult& result, int meshX, int meshY) {
    FieldUtil::InitializeField(result.f.u, meshX, meshY, 0);
    FieldUtil::InitializeField(result.f.v, meshX, meshY, 0);
    setInflowBoundaryCondition(result.f, meshX, meshY);
    FieldUtil::InitializeField(result.p, meshX, meshY, 0);
    FieldUtil::InitializeField(result.s, meshX, meshY, 0);
    FieldUtil::InitializeField(result.rot, meshX, meshY, 0);
    result.drag.x = 0.0;
    result.drag.y = 0.0;
}

void executeApplication() {
    std::filesystem::path configFile = std::filesystem::current_path() / "flow.h5";
    Hdf5Util file(configFile.string());
    int meshX, meshY, poissonIteration, interval, maxIterations;
    Value dx, dy, dt, reynolds, omega, epsilon, pRef;
    readConfig(file, meshX, meshY, dx, dy, dt, reynolds, omega, epsilon, pRef, 
        poissonIteration, maxIterations, interval);
    MeshRange2d range = {1, meshX - 3, 1, meshY - 3};

    AnalysisResult result;
    InitializeAnalysisResult(result, meshX, meshY);

    Object object;
    defineObject(object, meshX, meshY);

    NavierStokes2d solver(meshX, meshY, reynolds, dx, dy, dt,
        omega, epsilon, pRef, poissonIteration, range, result, object);
    for (int time = 1; time <= maxIterations; time++) {
        result = solver.calculate();

        if (time % interval != 0) {
            continue;
        }

        file.saveResult(result.rot, "rot", time);
        printf("dragX = %6.3f, dragY = %6.3f\n", result.drag.x, result.drag.y);
    }
}

int main() {
    try {
        executeApplication();
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}