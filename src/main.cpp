#include <iostream>
#include <iomanip>
#include <math.h>
#include "diff.h"

using namespace std;

// dx/dt = (3 - 2y)x
/* float eq1(float t, float x, float y) {
    return (3 - 2 * y) * x;
} */

//dy/dt = (2 - 2x)y
/* float eq2(float t, float x, float y) {
    return (2 * x - 2) * y;
} */

// d2y/dt2 = - g/l sin(y) g=9.81 l=1

float pendulum_eq(float t, float y, float dy) {
    return -(9.81 / 1) * sin(y);
}

int main() {
    /* FirstOrderODE* p = new FirstOrderODE(0, 15, 1);
    p->solve(&diff, 0.1);
    p->print_sol(); */

    /* FirstOrderODETwoSystem* p = new FirstOrderODETwoSystem(0, 10, 1, 1);
    p->solve(&eq1, &eq2, 0.001);
    p->print_sol();
    p->save_to_csv("lotka_volterra.csv"); */

    SecondOrderODE* p = new SecondOrderODE(0, 10, 3.054, 0);
    p->solve(&pendulum_eq, 0.001);
    p->print_sol();
    p->save_to_csv("pendulum_20cm.csv");

    return 0;
}
