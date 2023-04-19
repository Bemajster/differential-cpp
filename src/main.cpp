#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

float eq(float t, float y, float dy, float ddy) {
    return cos(y) + sin(t);
}

int main() {
    ThirdOrderODE* p = new ThirdOrderODE(0, 10, 1, 0, 0);
    p->solveRK4(&eq, 0.01);
    p->save_to_csv("differential.csv");

    return 0;
}
