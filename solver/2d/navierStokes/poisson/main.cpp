#include "Poisson2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>
#include <math.h>

void initializeSourceTerm(Field2d& s, int meshX, int meshY, double dx, double dy) {
    const double PI=3.14159;
    const double K_X = 2 * PI;
    const double K_Y = 2 * PI;
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            s[j][i] = -1 * (K_X * K_X + K_Y * K_Y) * sin(K_X * i * dx) * sin(K_Y * j * dy);
        }
    }
}

int main() {
    int meshX = 10;
    int meshY = meshX;
    double lx = 1.0;
    double ly = lx;
    double dx = lx / (meshX - 1);
    double dy = ly / (meshY - 1);
    double omega = 1.4;
    double epsilon = 1e-7;
    double pRef = 1.0;
    MeshRange2d range = {1, meshX - 2, 1, meshY - 2};
    Field2d s, p;
    FieldUtil::setSize(s, meshX, meshY);
    FieldUtil::setSize(p, meshX, meshY);
    initializeSourceTerm(s, meshX, meshY, dx, dy);

    try {
        Poisson2d poisson(meshX, meshY, dx, dy, omega, epsilon, pRef, range);
        const int interval = 1;
        const int maxIterations = 100;

        FlagField ip;
        FieldUtil::InitializeFlagField(ip, meshX, meshY, ObjectFlag::fluid);

        int iteration = poisson.calculate(p, s, ip, maxIterations);
        FieldUtil::display(p, iteration, interval);

        printf("Iteration: %d\n", iteration);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}