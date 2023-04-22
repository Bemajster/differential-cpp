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

float func_list(int n, float t, float x[]) {
    switch(n) {
        case 0:
            return 10 * (x[1] - x[0]);
            break;
        case 1:
            return x[0]*(28-x[2])-x[1];
            break;
        case 2:
            return x[0] * x[1] - 8 * x[2] / 3;
            break;
        default:
            return 0;
            break;
    }
}

int main() {
    /*FirstOrderODEFourSystem* p = new FirstOrderODEFourSystem(0, 5, 1, 0.5, 0, -1);

    p->solveRK4(&eq1, &eq2, &eq3, &eq4, 0.01);
    p->save_to_csv("four_RK4_system.csv");*/

    float sol0[] = {1, 0.5, 0};

    FirstOrderODESystem* s = new FirstOrderODESystem(0, 30, 3, sol0);
    s->solveRK2(&func_list, 0.001);
    s->print_sol();

    return 0;
}
