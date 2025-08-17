#include "Diffusion1d.h"
#include "FieldUtil.h"
#include <stdexcept>
#include <iostream>

void initializeField(Field1d& f, int meshX) {
    for (int i = 0; i < meshX; i++) {
        if (4 <= i && i <= 5) {
            f[i] = 1.0;
        }
        else {
            f[i] = 0.0;
        } 
    }
}

int main() {
    int meshX = 10;
    double constA = 1e-5;
    double deltaX = 0.001;
    double deltaT = 0.01;
    Field1d f;
    f.resize(meshX);
    initializeField(f, meshX);

    try {
        Diffusion1d diffusion(meshX, constA, deltaX, deltaT);
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