#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

float func_list(int n, float t, float x[]) {
    switch(n) {
        case 0:
            return (x[2] - 0.7) * x[0] - 3.5 * x[1];
            break;
        case 1:
            return 3.5 * x[0] + (x[2] - 0.7) * x[1];
            break;
        case 2:
            return 0.6 + 0.95 * x[2] - x[2] * x[2] * x[2] / 3 - (x[0] * x[0] + x[1] * x[1]) * (1 + 0.25 * x[2]) + 0.1 * x[2] * x[0] * x[0] * x[0];
            break;
        default:
            return 0;
            break;
    }
}

int main() {
    float sol0[] = {0.1, 0, 0};

    FirstOrderODESystem* s = new FirstOrderODESystem(0, 200, 3, sol0);
    s->solveRK2(&func_list, 0.01);
    s->save_to_csv("aizawa.csv", ";");

    return 0;
}
