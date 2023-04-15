#include <iostream>
#include <iomanip>
#include <math.h>
#include "diff.h"

using namespace std;

// dy/dx = diff(x, y)
float eq1(float t, float x, float y) {
    return (3 - 2 * y) * x;
}

float eq2(float t, float x, float y) {
    return (2 * x - 2) * y;
}


int main() {
    /* FirstOrderODE* p = new FirstOrderODE(0, 15, 1);
    p->solve(&diff, 0.1);
    p->print_sol(); */

    FirstOrderODETwoSystem* p = new FirstOrderODETwoSystem(0, 5, 1, 1);
    p->solve(&eq1, &eq2, 0.1);
    p->print_sol();
    p->save_to_csv("lotka_volterra.csv");

    return 0;
}
