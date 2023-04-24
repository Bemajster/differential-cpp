#include "ODE.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <array>

FirstOrderODE::FirstOrderODE(float p_a, float p_b, float p_y0) {
    a = p_a;
    b = p_b;
    y0 = p_y0;
}

FirstOrderODE::~FirstOrderODE() {
    delete sol;
}

void FirstOrderODE::solveRK2(float (*func)(float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol = new float[arr_size];

    sol[0] = y0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = step * func(a + i * step, sol[i]);
        float k2 = step * func(a + i * step + step, sol[i] + k1);

        sol[i+1] = sol[i] + (k1 + k2) / 2;
    }
}

void FirstOrderODE::solveRK4(float (*func)(float, float), float p_step) {
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

void FirstOrderODE::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol[i] << std::endl;
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

void FirstOrderODETwoSystem::solveRK2(float (*func_x)(float, float, float), float (*func_y)(float, float, float), float p_step) {
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

void FirstOrderODETwoSystem::solveRK4(float (*func_x)(float, float, float), float (*func_y)(float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol_x = new float[arr_size];
    sol_y = new float[arr_size];

    sol_x[0] = x0;
    sol_y[0] = y0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = func_x(a + i * step, sol_x[i], sol_y[i]);
        float l1 = func_y(a + i * step, sol_x[i], sol_y[i]);

        float k2 = func_x(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2);
        float l2 = func_y(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2);

        float k3 = func_x(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2);
        float l3 = func_y(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2);

        float k4 = func_x(a + i * step + step, sol_x[i] + step * k3, sol_y[i] + step * l3);
        float l4 = func_y(a + i * step + step, sol_x[i] + step * k3, sol_y[i] + step * l3);

        sol_x[i+1] = sol_x[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
        sol_y[i+1] = sol_y[i] + (l1 + 2 * l2 + 2 * l3 + l4) * step / 6;
    }
}

void FirstOrderODETwoSystem::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "t = " << a + i * step << " => x = " << sol_x[i] << " y = " << sol_y[i] << '\n';
    }
}

void FirstOrderODETwoSystem::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol_x[i] << separator << sol_y[i] << std::endl;
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

void SecondOrderODE::solveRK2(float (*func)(float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float *sol_dy;

    sol = new float[arr_size];
    sol_dy = new float[arr_size];

    sol[0] = y0;
    sol_dy[0] = dy0;

    for(int i = 0; i < arr_size; i++) {
        float k1 = step * sol_dy[i];
        float l1 = step * func(a + i * step, sol[i], sol_dy[i]);

        float k2 = step * (sol_dy[i] + l1);
        float l2 = step * func(a + i * step + step, sol[i] + k1, sol_dy[i] + l1);

        sol[i+1] = sol[i] + (k1 + k2) / 2;
        sol_dy[i+1] = sol_dy[i] + (l1 + l2) / 2;
    }

    delete sol_dy;
}

void SecondOrderODE::solveRK4(float (*func)(float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float *sol_dy;

    sol = new float[arr_size];
    sol_dy = new float[arr_size];

    sol[0] = y0;
    sol_dy[0] = dy0;

    for(int i = 0; i < arr_size; i++) {
        float k1 = sol_dy[i];
        float l1 = func(a + i * step, sol[i], sol_dy[i]);

        float k2 = sol_dy[i] + step * k1 / 2;
        float l2 = func(a + i * step + step / 2, sol[i] + step * k1 / 2, sol_dy[i] + step * l1 / 2);

        float k3 = sol_dy[i] + step * k2 / 2;
        float l3 = func(a + i * step + step / 2, sol[i] + step * k2 / 2, sol_dy[i] + step * l2 / 2);

        float k4 = sol_dy[i] + step * k3;
        float l4 = func(a + i * step + step, sol[i] + step * k3, sol_dy[i] + step * l3);

        sol[i+1] = sol[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
        sol_dy[i+1] = sol_dy[i] + (l1 + 2 * l2 + 2 * l3 + l4) * step / 6;
    }

    delete sol_dy;
}

void SecondOrderODE::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "x = " << a + i * step << " => y = " << sol[i] << '\n';
    }
}

void SecondOrderODE::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol[i] << std::endl;
    }

    output.close();
}

FirstOrderODEThreeSystem::FirstOrderODEThreeSystem(float p_a, float p_b, float p_x0, float p_y0, float p_z0) {
    a = p_a;
    b = p_b;
    x0 = p_x0;
    y0 = p_y0;
    z0 = p_z0;
}

FirstOrderODEThreeSystem::~FirstOrderODEThreeSystem() {
    delete sol_x;
    delete sol_y;
    delete sol_z;
}

void FirstOrderODEThreeSystem::solveRK2(float (*func_x)(float, float, float, float), float (*func_y)(float, float, float, float), float(*func_z)(float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol_x = new float[arr_size];
    sol_y = new float[arr_size];
    sol_z = new float[arr_size];

    sol_x[0] = x0;
    sol_y[0] = y0;
    sol_z[0] = z0;

    for(int i = 0; i < arr_size; i++) {
        float k1 = step * func_x(a + i * step, sol_x[i], sol_y[i], sol_z[i]);
        float l1 = step * func_y(a + i * step, sol_x[i], sol_y[i], sol_z[i]);
        float m1 = step * func_z(a + i * step, sol_x[i], sol_y[i], sol_z[i]);

        float k2 = step * func_x(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1);
        float l2 = step * func_y(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1);
        float m2 = step * func_z(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1);

        sol_x[i+1] = sol_x[i] + (k1 + k2) / 2;
        sol_y[i+1] = sol_y[i] + (l1 + l2) / 2;
        sol_z[i+1] = sol_z[i] + (m1 + m2) / 2;
    }
}

void FirstOrderODEThreeSystem::solveRK4(float (*func_x)(float, float, float, float), float (*func_y)(float, float, float, float), float(*func_z)(float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol_x = new float[arr_size];
    sol_y = new float[arr_size];
    sol_z = new float[arr_size];

    sol_x[0] = x0;
    sol_y[0] = y0;
    sol_z[0] = z0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = func_x(a + i * step, sol_x[i], sol_y[i], sol_z[i]);
        float l1 = func_y(a + i * step, sol_x[i], sol_y[i], sol_z[i]);
        float m1 = func_z(a + i * step, sol_x[i], sol_y[i], sol_z[i]);

        float k2 = func_x(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2);
        float l2 = func_y(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2);
        float m2 = func_z(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2);
    
        float k3 = func_x(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2);
        float l3 = func_y(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2);
        float m3 = func_z(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2);

        float k4 = func_x(a + i * step + step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3);
        float l4 = func_y(a + i * step + step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3);
        float m4 = func_z(a + i * step + step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3);
    
        sol_x[i+1] = sol_x[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
        sol_y[i+1] = sol_y[i] + (l1 + 2 * l2 + 2 * l3 + l4) * step / 6;
        sol_z[i+1] = sol_z[i] + (m1 + 2 * m2 + 2 * m3 + m4) * step / 6;
    }
}

void FirstOrderODEThreeSystem::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "t = " << a + i * step << " => x = " << sol_x[i] << " y = " << sol_y[i] << " z = " << sol_z[i] << '\n';
    }
}

void FirstOrderODEThreeSystem::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol_x[i] << separator << sol_y[i] << separator << sol_z[i] << std::endl;
    }

    output.close();
}

ThirdOrderODE::ThirdOrderODE(float p_a, float p_b, float p_y0, float p_dy0, float p_ddy0) {
    a = p_a;
    b = p_b;
    y0 = p_y0;
    dy0 = p_dy0;
    ddy0 = p_ddy0;
}

ThirdOrderODE::~ThirdOrderODE() {
    delete sol;
}

void ThirdOrderODE::solveRK2(float (*func)(float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float *sol_dy, *sol_ddy;

    sol = new float[arr_size];
    sol_dy = new float[arr_size];
    sol_ddy = new float[arr_size];

    sol[0] = y0;
    sol_dy[0] = dy0;
    sol_ddy[0] = ddy0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = step * sol_dy[i];
        float l1 = step * sol_ddy[i];
        float m1 = step * func(a + i * step, sol[i], sol_dy[i], sol_ddy[i]);

        float k2 = step * (sol_dy[i] + l1);
        float l2 = step * (sol_ddy[i] + m1);
        float m2 = step * func(a + i * step + step, sol[i] + k1, sol_dy[i] + l1, sol_ddy[i] + m1);

        sol[i+1] = sol[i] + (k1 + k2) / 2;
        sol_dy[i+1] = sol_dy[i] + (l1 + l2) / 2;
        sol_ddy[i+1] = sol_ddy[i] + (m1 + m2) / 2;
    }

    delete sol_dy;
    delete sol_ddy;
}

void ThirdOrderODE::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "x = " << a + i * step << " => y = " << sol[i] << '\n';
    }
}

void ThirdOrderODE::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol[i] << std::endl;
    }

    output.close();
}

void ThirdOrderODE::solveRK4(float (*func)(float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float *sol_dy, *sol_ddy;

    sol = new float[arr_size];
    sol_dy = new float[arr_size];
    sol_ddy = new float[arr_size];

    sol[0] = y0;
    sol_dy[0] = dy0;
    sol_ddy[0] = ddy0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = sol_dy[i];
        float l1 = sol_ddy[i];
        float m1 = func(a + i * step, sol[i], sol_dy[i], sol_ddy[i]);

        float k2 = sol_dy[i] + step * l1 / 2;
        float l2 = sol_ddy[i] + step * m1 / 2;
        float m2 = func(a + i * step + step / 2, sol[i] + step * k1 / 2, sol_dy[i] + step * l1 / 2, sol_ddy[i] + step * m1 / 2);

        float k3 = sol_dy[i] + step * l2 / 2;
        float l3 = sol_ddy[i] + step * m2 / 2;
        float m3 = func(a + i * step + step / 2, sol[i] + step * k2 / 2, sol_dy[i] + step * l2 / 2, sol_ddy[i] + step * m2 / 2);

        float k4 = sol_dy[i] + step * k3;
        float l4 = sol_ddy[i] + step * l3;
        float m4 = func(a + i * step + step, sol[i] + step * k3, sol_dy[i] + step * l3, sol_ddy[i] + step * m3);

        sol[i+1] = sol[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
        sol_dy[i+1] = sol_dy[i] + (l1 + 2 * l2 + 2 * l3 + l4) * step / 6;
        sol_ddy[i+1] = sol_ddy[i] + (m1 + 2 * m2 + 2 * m3 + m4) * step / 6;
    }

    delete sol_dy;
    delete sol_ddy;
}

FirstOrderODEFourSystem::FirstOrderODEFourSystem(float p_a, float p_b, float p_x0, float p_y0, float p_z0, float p_w0) {
    a = p_a;
    b = p_b;
    x0 = p_x0;
    y0 = p_y0;
    z0 = p_z0;
    w0 = p_w0;
}

FirstOrderODEFourSystem::~FirstOrderODEFourSystem() {
    delete sol_x;
    delete sol_y;
    delete sol_z;
    delete sol_w;
}

void FirstOrderODEFourSystem::solveRK2(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float(*func_z)(float, float, float, float, float), float (*func_w)(float, float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol_x = new float[arr_size];
    sol_y = new float[arr_size];
    sol_z = new float[arr_size];
    sol_w = new float[arr_size];

    sol_x[0] = x0;
    sol_y[0] = y0;
    sol_z[0] = z0;
    sol_w[0] = w0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = step * func_x(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);
        float l1 = step * func_y(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);
        float m1 = step * func_z(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);
        float n1 = step * func_w(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);

        float k2 = step * func_x(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1, sol_w[i] + n1);
        float l2 = step * func_y(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1, sol_w[i] + n1);
        float m2 = step * func_z(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1, sol_w[i] + n1);
        float n2 = step * func_w(a + i * step + step, sol_x[i] + k1, sol_y[i] + l1, sol_z[i] + m1, sol_w[i] + n1);

        sol_x[i+1] = sol_x[i] + (k1 + k2) / 2;
        sol_y[i+1] = sol_y[i] + (l1 + l2) / 2;
        sol_z[i+1] = sol_z[i] + (m1 + m2) / 2;
        sol_w[i+1] = sol_w[i] + (n1 + n2) / 2;
    }
}

void FirstOrderODEFourSystem::solveRK4(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float(*func_z)(float, float, float, float, float), float (*func_w)(float, float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol_x = new float[arr_size];
    sol_y = new float[arr_size];
    sol_z = new float[arr_size];
    sol_w = new float[arr_size];

    sol_x[0] = x0;
    sol_y[0] = y0;
    sol_z[0] = z0;
    sol_w[0] = w0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = func_x(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);
        float l1 = func_y(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);
        float m1 = func_z(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);
        float n1 = func_w(a + i * step, sol_x[i], sol_y[i], sol_z[i], sol_w[i]);

        float k2 = func_x(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2, sol_w[i] + step * n1 / 2);
        float l2 = func_y(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2, sol_w[i] + step * n1 / 2);
        float m2 = func_z(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2, sol_w[i] + step * n1 / 2);
        float n2 = func_w(a + i * step + step / 2, sol_x[i] + step * k1 / 2, sol_y[i] + step * l1 / 2, sol_z[i] + step * m1 / 2, sol_w[i] + step * n1 / 2);
    
        float k3 = func_x(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2, sol_w[i] + step * n2 / 2);
        float l3 = func_y(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2, sol_w[i] + step * n2 / 2);
        float m3 = func_z(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2, sol_w[i] + step * n2 / 2);
        float n3 = func_w(a + i * step + step / 2, sol_x[i] + step * k2 / 2, sol_y[i] + step * l2 / 2, sol_z[i] + step * m2 / 2, sol_w[i] + step * n2 / 2);
    
        float k4 = func_x(a + i * step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3, sol_w[i] + step * n3);
        float l4 = func_y(a + i * step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3, sol_w[i] + step * n3);
        float m4 = func_z(a + i * step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3, sol_w[i] + step * n3);
        float n4 = func_w(a + i * step, sol_x[i] + step * k3, sol_y[i] + step * l3, sol_z[i] + step * m3, sol_w[i] + step * n3);

        sol_x[i+1] = sol_x[i] + (k1 + 2 * k2 + 2 * k3 + k4) * step / 6;
        sol_y[i+1] = sol_y[i] + (l1 + 2 * l2 + 2 * l3 + l4) * step / 6;
        sol_z[i+1] = sol_z[i] + (m1 + 2 * m2 + 2 * m3 + m4) * step / 6;
        sol_w[i+1] = sol_w[i] + (n1 + 2 * n2 + 2 * n3 + n4) * step / 6;
    }
}

void FirstOrderODEFourSystem::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "t = " << a + i * step << " => x = " << sol_x[i] << " y = " << sol_y[i] << " z = " << sol_z[i] << " w = " << sol_w[i] << '\n';
    }
}

void FirstOrderODEFourSystem::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol_x[i] << separator << sol_y[i] << separator << sol_z[i] << separator << sol_w[i] << std::endl;
    }

    output.close();
}

SecondOrderODETwoSystem::SecondOrderODETwoSystem(float p_a, float p_b, float p_x0, float p_dx0, float p_y0, float p_dy0) {
    a = p_a;
    b = p_b;
    x0 = p_x0;
    dx0 = p_dx0;
    y0 = p_y0;
    dy0 = p_dy0;
}

SecondOrderODETwoSystem::~SecondOrderODETwoSystem() {
    delete sol_x;
    delete sol_y;
}

void SecondOrderODETwoSystem::solveRK2(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float *sol_dx, *sol_dy;

    sol_x = new float[arr_size];
    sol_dx = new float[arr_size];
    sol_y = new float[arr_size];
    sol_dy = new float[arr_size];

    sol_x[0] = x0;
    sol_dx[0] = dx0;
    sol_y[0] = y0;
    sol_dy[0] = dy0;

    for(int i = 0; i <= arr_size - 1; i++) {
        float k1 = step * sol_dx[i];
        float l1 = step * sol_dy[i];
        float m1 = step * func_x(a + i * step, sol_x[i], sol_dx[i], sol_y[i], sol_dy[i]);
        float n1 = step * func_y(a + i * step, sol_x[i], sol_dx[i], sol_y[i], sol_dy[i]);

        float k2 = step * (sol_dx[i] + l1);
        float l2 = step * (sol_dy[i] + n1);
        float m2 = step * func_x(a + i * step + step, sol_x[i] + k1, sol_dx[i] + l1, sol_y[i] + m1, sol_dy[i] + n1);
        float n2 = step * func_y(a + i * step + step, sol_x[i] + k1, sol_dx[i] + l1, sol_y[i] + m1, sol_dy[i] + n1);

        sol_x[i+1] = sol_x[i] + (k1 + k2) / 2;
        sol_y[i+1] = sol_y[i] + (l1 + l2) / 2;
        sol_dx[i+1] = sol_dx[i] + (m1 + m2) / 2;
        sol_dy[i+1] = sol_dy[i] + (n1 + n2) / 2; 
    }

    delete sol_dx;
    delete sol_dy;
}

void SecondOrderODETwoSystem::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "t = " << a + i * step << " => x = " << sol_x[i] << " y = " << sol_y[i] << '\n';
    }
}

void SecondOrderODETwoSystem::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol_x[i] << separator << sol_y[i] << std::endl;
    }

    output.close();
}

FirstOrderODESystem::FirstOrderODESystem(float p_a, float p_b, int p_order, float p_sol0[]) {
    a = p_a;
    b = p_b;

    order = p_order;

    sol0 = new float[order];

    sol0 = p_sol0;
}

FirstOrderODESystem::~FirstOrderODESystem() {
    delete sol0;
    
    for(int i = 0; i < order; i++) {
        delete sol[i];
    }

    delete sol;
}

void FirstOrderODESystem::solveRK2(float (*func_list)(int, float, float[]), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    sol = new float*[order];

    for(int i = 0; i < order; i++) {
        sol[i] = new float[arr_size];
        sol[i][0] = sol0[i];
    }

    for(int i = 0; i < arr_size; i++) {
        float k1n[order], k2n[order];

        for(int j = 0; j < order; j++) {
            float args[order];

            for(int k = 0; k < order; k++) {
                args[k] = sol[k][i];
            }

            k1n[j] = step * func_list(j, a + i * step, args);
        }

        for(int j = 0; j < order; j++) {
            float args[order];

            for(int k = 0; k < order; k++) {
                args[k] = sol[k][i] + k1n[k];
            }

            k2n[j] = step * func_list(j, a + i * step + step, args);
        }

        for(int j = 0; j < order; j++) {
            sol[j][i + 1] = sol[j][i] + (k1n[j] + k2n[j]) / 2;
        }
    }
}

void FirstOrderODESystem::solveRK4(float (*func_list)(int, float, float[]), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    for(int i = 0; i < order; i++) {
        sol[i] = new float[arr_size];
        sol[i][0] = sol0[i];
    }

    for(int i = 0; i < arr_size; i++) {
        float k1n[order], k2n[order], k3n[order], k4n[order];

        for(int j = 0; j < order; j++) {
            float args[order];

            for(int k = 0; k < order; k++) {
                args[k] = sol[k][i];
            }

            k1n[j] = func_list(j, a + i * step, args);
        }

        for(int j = 0; j < order; j++) {
            float args[order];

            for(int k = 0; k < order; k++) {
                args[k] = sol[k][i] + step * k1n[k] / 2;
            }

            k2n[j] = func_list(j, a + i * step + step / 2, args);
        }

        for(int j = 0; j < order; j++) {
            float args[order];

            for(int k = 0; k < order; k++) {
                args[k] = sol[k][i] + step * k2n[k] / 2;
            }

            k3n[j] = func_list(j, a + i * step + step / 2, args);
        }

        for(int j = 0; j < order; j++) {
            float args[order];

            for(int k = 0; k < order; k++) {
                args[k] = sol[k][i] + step * k3n[k];
            }

            k4n[j] = func_list(j, a + i * step, args);
        }

        for(int j = 0; j < order; j++) {
            sol[j][i + 1] = sol[j][i] + (k1n[j] + 2 * k2n[j] + 2 * k3n[j] + k4n[j]) * step / 6;
        }
    }
}

void FirstOrderODESystem::print_sol() {
    for(int i = 0; i <= ceil((b - a) / step); i++) {
        std::cout << "t = " << a + i * step << " =>";
        for(int j = 0; j < order; j++) {
            std::cout << " x" << j << " = " << sol[j][i]; 
        }
        std::cout << '\n';
    }
}

void FirstOrderODESystem::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step;
        for(int j = 0; j < order; j++) {
            output << separator << sol[j][i];
        }
        output << std::endl;
    }

    output.close();
}

int FirstOrderODESystem::get_order() {
    return order;
}

ODE::ODE(float p_a, float p_b, int p_order, float p_sol0[]) {
    a = p_a; 
    b = p_b;
    order = p_order;

    sol0 = new float[order];

    sol0 = p_sol0;
}

ODE::~ODE() {
    delete sol0;
    delete sol;
}

void ODE::solveRK2(float (*func)(float, float[]), float p_step) {
    step = p_step;

    int arr_size = ceil((b - a) / step) + 1;

    float **sol_d;
    sol_d = new float*[order];

    for(int i = 0; i < order; i++) {
        sol_d[i] = new float[arr_size];
        sol_d[i][0] = sol0[i];
    }

    for(int i = 0; i < arr_size; i++) {
        float k1n[order], k2n[order], args1[order], args2[order];

        for(int j = 0; j < order - 1; j++) {
            k1n[j] = step * sol_d[j+1][i];
        }

        for(int j = 0; j < order; j++) {
            args1[j] = sol_d[j][i];
        }

        k1n[order - 1] = step * func(a + i * step, args1);

        for(int j = 0; j < order - 1; j++) {
            k2n[j] = step * (sol_d[j+1][i] + k1n[j]);
        }

        for(int j = 0; j < order; j++) {
            args2[j] = sol_d[j][i] + k1n[j];
        }

        k2n[order - 1] = step * func(a + i * step + step, args2);

        for(int j = 0; j < order; j++) {
            sol_d[j][i + 1] = sol_d[j][i] + (k1n[j] + k2n[j]) / 2;
        }
    }

    sol = new float[arr_size];
    sol = sol_d[0];

    delete sol_d;
}

void ODE::save_to_csv(std::string dir, std::string separator) {
    std::fstream output;
    output.open(dir, std::ios::out);

    for(int i = 0; i <= ceil((b - a) / step); i++) {
        output << a + i * step << separator << sol[i] << std::endl;
    }

    output.close();
}