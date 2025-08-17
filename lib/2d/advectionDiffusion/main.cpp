#include "AdvectionDiffusion2d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void initializeField(Field2d& f, int meshX, int meshY) {
    for (int j = 0; j < meshY; j++) {
        for (int i = 0; i < meshX; i++) {
            if (4 <= i && i <= 5 &&
                4 <= j && j <= 5) {
                f[j][i] = 1.0;
            }
            else {
                f[j][i] = 0.0;
            }
        }
    }
}

int main() {
    int meshX = 10;
    int meshY = meshX;
    double constA = 1e-5;
    double constU = 0.02;
    double constV = 0.00;
    double deltaX = 0.001;
    double deltaY = deltaX;
    double deltaT = 0.01; 
    Field2d f;
    FieldUtil::setSize(f, meshX, meshY);
    initializeField(f, meshX, meshY);

    try {
        AdvectionDiffusion2d diffusion(meshX, meshY, constA, constU, constV, deltaX, deltaY, deltaT);
        const int interval = 1;
        const int maxIterations = 10;
        for (int time = 0; time < maxIterations; time++) {
            FieldUtil::display(f, time, interval);
            f = diffusion.calculate(f);
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}