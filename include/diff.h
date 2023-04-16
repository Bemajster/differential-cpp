#pragma once
#include <iostream>

class FirstOrderODE {
public:
    FirstOrderODE(float a, float b, float y0);
    ~FirstOrderODE();

    void solve(float (*func)(float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

private:
    float a, b, y0, step, *sol;
};

class SecondOrderODE {
public:
    SecondOrderODE(float a, float b, float y0, float dy0);
    ~SecondOrderODE();

    void solve(float (*func)(float, float, float), float step); // f(x, y, dy)
    void print_sol();
    void save_to_csv(std::string dir);

private:
    float a, b, y0, dy0, step, *sol;
};

class FirstOrderODETwoSystem {
public:
    FirstOrderODETwoSystem(float a, float b, float x0, float y0);
    ~FirstOrderODETwoSystem();

    void solve(float (*func_x)(float, float, float), float (*func_y)(float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

private:
    float a, b, x0, y0, step, *sol_x, *sol_y;
};

class FirstOrderODEThreeSystem {
public:
    FirstOrderODEThreeSystem(float a, float b, float x0, float y0, float z0);
    ~FirstOrderODEThreeSystem();

    void solveRK2(float (*func_x)(float, float, float, float), float (*func_y)(float, float, float, float), float(*func_z)(float, float, float, float), float step);
    void solveRK4(float (*func_x)(float, float, float, float), float (*func_y)(float, float, float, float), float(*func_z)(float, float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

private:
    float a, b, x0, y0, z0, step, *sol_x, *sol_y, *sol_z;
};