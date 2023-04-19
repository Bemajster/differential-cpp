#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

float eq1(float t, float x, float dx, float y, float dy) {
    return -y;
}

float eq2(float t, float x, float dx, float y, float dy) {
    return -x;
}

int main() {
    SecondOrderODETwoSystem* p = new SecondOrderODETwoSystem(0, 10, 1, -1, 0, 0);
    p->solveRK2(&eq1, &eq2, 0.1);
    p->print_sol();
    p->save_to_csv("second_order_system.csv");

    return 0;
}
