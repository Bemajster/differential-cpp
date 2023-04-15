#pragma once
#include <iostream>

class FirstOrderDiff {
public:
    FirstOrderDiff(float a, float b, float y0);
    ~FirstOrderDiff();

    void solve(float (*func)(float, float), float step);
    void print_sol();
    void save_to_csv(std::string dir);

private:
    float a, b, y0, step;
    float *sol;
};
