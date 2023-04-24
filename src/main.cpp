#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

float eq(float t, float x[]) {
    return -9.81 * sin(x[0]);
}

int main() {
    float sol0[] = {3.05, 0};

    ODE* s = new ODE(0, 10, 2, sol0);
    s->solveRK2(&eq, 0.0001);
    s->save_to_csv("pendulum.csv", ";");

    return 0;
}
