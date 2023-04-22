#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

/*float eq1(float t, float x, float y, float z, float w) {
    return (5 - x) * y;
}

float eq2(float t, float x, float y, float z, float w) {
    return - z * x;
}

float eq3(float t, float x, float y, float z, float w) {
    return (8 + w) - y;
}

float eq4(float t, float x, float y, float z, float w) {
    return y-z;
}*/

float eq(float t, float y) {
    return y;
}

int main() {
    /*FirstOrderODEFourSystem* p = new FirstOrderODEFourSystem(0, 5, 1, 0.5, 0, -1);

    p->solveRK4(&eq1, &eq2, &eq3, &eq4, 0.01);
    p->save_to_csv("four_RK4_system.csv");*/

    FirstOrderODE* p = new FirstOrderODE(0, 5, 1);
    p->solveRK2(&eq, 0.01);
    p->save_to_csv("exponential.csv");

    return 0;
}
