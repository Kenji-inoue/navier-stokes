#include <gtest/gtest.h>
#include "Poisson2d.h"
#include "FieldUtil.h"

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

TEST(poissonTest, simulation) {
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

    Poisson2d poisson(meshX, meshY, dx, dy, omega, epsilon, pRef, range);
    const int interval = 1;
    const int maxIterations = 100;

    FlagField ip;
    FieldUtil::InitializeFlagField(ip, meshX, meshY, ObjectFlag::fluid);
    int iteration = poisson.calculate(p, s, ip, maxIterations);

    std::vector<std::vector<double>> expected = {
    { 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000 },
    { 0.000,  0.430,  0.659,  0.580,  0.229, -0.229, -0.580, -0.659, -0.430,  0.000 },
    { 0.000,  0.659,  1.010,  0.888,  0.351, -0.351, -0.888, -1.010, -0.659,  0.000 },
    { 0.000,  0.580,  0.888,  0.781,  0.309, -0.309, -0.781, -0.888, -0.580,  0.000 },
    { 0.000,  0.229,  0.351,  0.309,  0.122, -0.122, -0.309, -0.351, -0.229,  0.000 },
    { 0.000, -0.229, -0.351, -0.309, -0.122,  0.122,  0.309,  0.351,  0.229,  0.000 },
    { 0.000, -0.580, -0.888, -0.781, -0.309,  0.309,  0.781,  0.888,  0.580,  0.000 },
    { 0.000, -0.659, -1.010, -0.888, -0.351,  0.351,  0.888,  1.010,  0.659,  0.000 },
    { 0.000, -0.430, -0.659, -0.580, -0.229,  0.229,  0.580,  0.659,  0.430,  0.000 },
    { 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000 }
};

   for (int j = 0; j < meshY; j++) {
       for (int i = 0; i < meshX; i++) {
            printf("p[%d][%d] = %f\n", j, i, p[j][i]);
           EXPECT_NEAR(p[j][i], expected[j][i], 1e-3);
       }
   }

   EXPECT_EQ(iteration, 42);
}