#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

float func_list(int n, float t, float x[]) {
    switch(n) {
        case 0:
            return -x[0] + x[1] + x[1] * x[2];
            break;
        case 1:
            return -x[0] - x[1] + 0.4 * x[0] * x[2];
            break;
        case 2:
            return x[2] - 0.3 * x[0] * x[1];
            break;
        default:
            return 0;
            break;
    }
}

float *t, *sols;

int main() {
    float sol0[] = {1, -1, 1};

    FirstOrderODESystem* s = new FirstOrderODESystem(0, 200, 3, sol0);
    s->solveRK4(&func_list, 0.001);
    s->get_sol(t, sols);
    s->save_to_csv("Sakarya.csv", ",");
    cout << t[5] << endl;

    return 0;
}
