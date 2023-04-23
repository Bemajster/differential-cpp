#include <iostream>
#include <math.h>
#include "ODE.h"

using namespace std;

float eq1(float t, float x, float dx, float y, float dy) {
    return -y;
}

float eq2(float t, float x, float dx, float y, float dy) {
    return x;
}

int main() {
    SecondOrderODETwoSystem* s = new SecondOrderODETwoSystem(0, 10, 1, 0, 0, 0);
    s->solveRK2(&eq1, &eq2, 0.001);
    s->save_to_csv("second_order_sys_RK2.csv", ";");

    SecondOrderODETwoSystem* p = new SecondOrderODETwoSystem(0, 10, 1, 0, 0, 0);
    p->solveRK4(&eq1, &eq2, 0.001);
    p->save_to_csv("second_order_sys_RK4.csv", ";");

    return 0;
}
