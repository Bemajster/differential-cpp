#include <iostream>
#include <math.h>
#include "diff.h"

using namespace std;

float sigma = 10, r = 28, b = 8 / 3;

float eq1(float t, float x, float y, float z) {
    return sigma * y - sigma * x;
}

float eq2(float t, float x, float y, float z) {
    return -x * z + r * x - y;
}

float eq3(float t, float x, float y, float z) {
    return x * y - b * z;
}

float *sol_x, *sol_y, *sol_z, *t;

int main() {
    FirstOrderODEThreeSystem* p = new FirstOrderODEThreeSystem(0, 3000, 0, 0.5, 1);
    p->solveRK4(&eq1, &eq2, &eq3, 0.001);
    p->get_sol(t, sol_x, sol_y, sol_z);

    return 0;
}
