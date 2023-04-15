#include "diff.h"
#include <math.h>
#include <iostream>
#include <fstream>

FirstOrderDiff::FirstOrderDiff(float p_a, float p_b, float p_y0) {
    a = p_a;
    b = p_b;
    y0 = p_y0;
}

FirstOrderDiff::~FirstOrderDiff() {
    delete sol;
}

void FirstOrderDiff::solve(float (*func)(float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol = new float[arr_size];

    sol[0] = y0;

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        float k1 = func(a + i * step, sol[i]);
        float k2 = func(a + i * step + step / 2, sol[i] + step * k1 / 2);
        float k3 = func(a + i * step + step / 2, sol[i] + step * k2 / 2);
        float k4 = func(a + i * step + step, sol[i] + step * k3);

        sol[i+1] = sol[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
    }
}

void FirstOrderDiff::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "x = " << a + i * step << " => y = " << sol[i] << '\n';
    }
}

void FirstOrderDiff::save_to_csv(std::string dir) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << ";" << sol[i] << std::endl;
    }

    output.close();
}
