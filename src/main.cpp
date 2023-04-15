#include <iostream>
#include <math.h>
#include "diff.h"

using namespace std;

// dy/dx = diff(x, y)
float diff(float x, float y) {
    return x;
}

int main() {
    FirstOrderDiff* p = new FirstOrderDiff(0, 10, 0);
    p->solve(&diff, 5);
    p->print_sol();

    return 0;
}
