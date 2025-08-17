#include "Burgers2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void initializeVelocity(Velocity2d& f, int meshX, int meshY, double dx, double dy) {
    const double PI=3.14159;
    const double K = 2 * PI;
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            const auto x = i * dx;
            const auto y = j * dy;
            f.u[j][i] = -1 * cos(K * x) * sin(K * y);
            f.v[j][i] = sin(K * x) * cos(K * y);
        }
    }
}

int main() {
    int meshX = 10;
    int meshY = meshX;
    double reynolds = 200;
    double lx = 1.0;
    double ly = lx;
    double dx = lx / (meshX - 1);
    double dy = ly / (meshY - 1);
    double deltaT = 0.02; 
    Velocity2d f;
    FieldUtil::setSize(f.u, meshX, meshY);
    FieldUtil::setSize(f.v, meshX, meshY);
    initializeVelocity(f, meshX, meshY, dx, dy);

    try {
        Burgers2d burgers(meshX, meshY, reynolds, dx, dy, deltaT, f);
        const int interval = 1;
        const int maxIterations = 10;
        for (int time = 0; time < maxIterations; time++) {
            FieldUtil::display(f.u, time, interval);
            FieldUtil::display(f.v, time, interval);
            f = burgers.calculate();
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}