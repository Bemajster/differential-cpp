#include "diff.h"
#include <math.h>
#include <iostream>
#include <fstream>

FirstOrderODE::FirstOrderODE(float p_a, float p_b, float p_y0) {
    a = p_a;
    b = p_b;
    y0 = p_y0;
}

FirstOrderODE::~FirstOrderODE() {
    delete sol;
}

void FirstOrderODE::solve(float (*func)(float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol = new float[arr_size];

    sol[0] = y0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = func(a + i * step, sol[i]);
        float k2 = func(a + i * step + step / 2, sol[i] + step * k1 / 2);
        float k3 = func(a + i * step + step / 2, sol[i] + step * k2 / 2);
        float k4 = func(a + i * step + step, sol[i] + step * k3);

        sol[i+1] = sol[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
    }
}

void FirstOrderODE::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "x = " << a + i * step << " => y = " << sol[i] << '\n';
    }
}

void FirstOrderODE::save_to_csv(std::string dir) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << ";" << sol[i] << std::endl;
    }

    output.close();
}

FirstOrderODETwoSystem::FirstOrderODETwoSystem(float p_a, float p_b, float p_x0, float p_y0) {
    a = p_a;
    b = p_b;
    x0 = p_x0;
    y0 = p_y0;
}

FirstOrderODETwoSystem::~FirstOrderODETwoSystem() {
    delete sol_x;
    delete sol_y;
}

void FirstOrderODETwoSystem::solve(float (*func_x)(float, float, float), float (*func_y)(float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol_x = new float[arr_size];
    sol_y = new float[arr_size];

    sol_x[0] = x0;
    sol_y[0] = y0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = step * func_x(a + i * step, sol_x[i], sol_y[i]);
        float l1 = step * func_y(a + i * step, sol_x[i], sol_y[i]);

        float k2 = step * func_x(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1);
        float l2 = step * func_y(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1);

        sol_x[i+1] = sol_x[i] + (k1 + k2) / 2;
        sol_y[i+1] = sol_y[i] + (l1 + l2) / 2;
    }
}

void FirstOrderODETwoSystem::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "t = " << a + i * step << " => x = " << sol_x[i] << " y = " << sol_y[i] << '\n';
    }
}

void FirstOrderODETwoSystem::save_to_csv(std::string dir) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << ";" << sol_x[i] << ";" << sol_y[i] << std::endl;
    }

    output.close();
}

SecondOrderODE::SecondOrderODE(float p_a, float p_b, float p_y0, float p_dy0) {
    a = p_a;
    b = p_b;
    y0 = p_y0;
    dy0 = p_dy0;
}

SecondOrderODE::~SecondOrderODE() {
    delete sol;
}

void SecondOrderODE::solve(float (*func)(float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float *sol_dy;

    sol = new float[arr_size];
    sol_dy = new float[arr_size];

    sol[0] = y0;
    sol_dy[0] = dy0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = step * sol_dy[i];
        float l1 = step * func(a + i * step, sol[i], sol_dy[i]);

        float k2 = step * (sol_dy[i] + l1);
        float l2 = step * func(a + i * step + step, sol[i] + k1, sol_dy[i] + l1);

        sol[i+1] = sol[i] + (k1 + k2) / 2;
        sol_dy[i+1] = sol_dy[i] + (l1 + l2) / 2;
    }

    delete sol_dy;
}

void SecondOrderODE::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "x = " << a + i * step << " => y = " << sol[i] << '\n';
    }
}

void SecondOrderODE::save_to_csv(std::string dir) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << ";" << sol[i] << std::endl;
    }

    output.close();
}