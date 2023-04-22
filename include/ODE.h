#pragma once
#include <iostream>

class FirstOrderODE {
public:
    FirstOrderODE(float a, float b, float y0);
    ~FirstOrderODE();

    void solveRK2(float (*func)(float, float), float step);
    void solveRK4(float (*func)(float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float *t, float *y);

private:
    float a, b, y0, step, *sol;
};

class SecondOrderODE {
public:
    SecondOrderODE(float a, float b, float y0, float dy0);
    ~SecondOrderODE();

    void solveRK2(float (*func)(float, float, float), float step); // f(x, y, dy)
    void solveRK4(float (*func)(float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float *t, float *y);

private:
    float a, b, y0, dy0, step, *sol;
};

class ThirdOrderODE {
public:
    ThirdOrderODE(float a, float b, float y0, float dy0, float ddy0);
    ~ThirdOrderODE();

    void solveRK2(float (*func)(float, float, float, float), float step);
    void solveRK4(float (*func)(float, float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float *t, float *y);

private:
    float a, b, y0, dy0, ddy0, step, *sol;
};

class FirstOrderODETwoSystem {
public:
    FirstOrderODETwoSystem(float a, float b, float x0, float y0);
    ~FirstOrderODETwoSystem();

    void solveRK2(float (*func_x)(float, float, float), float (*func_y)(float, float, float), float step);
    void solveRK4(float (*func_x)(float, float, float), float (*func_y)(float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float *t, float *x, float *y);

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

    void get_sol(float *t, float *x, float *y, float *z);

private:
    float a, b, x0, y0, z0, step, *sol_x, *sol_y, *sol_z;
};

class FirstOrderODEFourSystem {
public:
    FirstOrderODEFourSystem(float a, float b, float x0, float y0, float z0, float w0);
    ~FirstOrderODEFourSystem();

    void solveRK2(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float(*func_z)(float, float, float, float, float), float (*func_w)(float, float, float, float, float), float step);
    void solveRK4(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float(*func_z)(float, float, float, float, float), float (*func_w)(float, float, float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float *t, float *x, float *y, float *z, float *w);

private:
    float a, b, x0, y0, z0, w0, step, *sol_x, *sol_y, *sol_z, *sol_w;
};

class SecondOrderODETwoSystem {
public:
    SecondOrderODETwoSystem(float a, float b, float x0, float dx0, float y0, float dy0);
    ~SecondOrderODETwoSystem();

    void solveRK2(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float step);
    void solveRK4(float (*func_x)(float, float, float, float, float), float (*func_y)(float, float, float, float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float *t, float *x, float *y);

private:
    float a, b, x0, dx0, y0, dy0, step, *sol_x, *sol_y;
};

class FirstOrderODESystem {
public:
    FirstOrderODESystem(float a, float b, float sol0[]);
    ~FirstOrderODESystem();

    void solveRK2(float (*func_list)(int, float[]), float step);
    void solveRK4(float (*func_list)(int, float[]), float step);
    void print_sol();
    void save_to_csv(std::string dir);

    void get_sol(float sols[]);

private:
    float a, b, step, **sol, *sol0;
    int order;
};