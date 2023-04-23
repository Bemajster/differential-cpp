#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;


float func_list(int n, float t, float x[]) {
    switch(n) {
        case 0:
            return 40 * (x[1] - x[0]) + 0.16 * x[0] * x[2];
            break;
        case 1:
            return 55 * x[0] + 20 * x[1] - x[0] * x[2];
            break;
        case 2:
            return 1.833 * x[2] + x[0] * x[1] - 0.65 * x[0] * x[0];
            break;
        default:
            return 0;
            break;
    }
}

float *t, *sols;

int main() {
    float sol0[] = {0.349, 0, -0.16};

    FirstOrderODESystem* s = new FirstOrderODESystem(0, 200, 3, sol0);
    s->solveRK4(&func_list, 0.001);
    s->save_to_csv("dequan-li.csv", ",");

    return 0;
}
