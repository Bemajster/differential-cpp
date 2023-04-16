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

int main() {
    FirstOrderODEThreeSystem* p = new FirstOrderODEThreeSystem(0, 50, 0, 0.5, 1);
    p->solve(&eq1, &eq2, &eq3, 0.0001);
    p->save_to_csv("lorenz_attractor.csv");

    return 0;
}
