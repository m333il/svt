#include <cstdio>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <fstream>
#include "inmost.h"

using namespace INMOST;

double f(double x, double y) {
    return sin(10 * x) * sin(10 * y);
}

double u(double x, double y, double dx = 1.0, double dy = 1.0) {
    return sin(10 * x) * sin(10 * y) / ((dx + dy) * 10 * 10);
}

double g(double x, double y) {
    return sin(10 * x) * sin(10 * y) / ((2) * 10 * 10);
}

int main(int argc, char **argv) {
    Solver::Initialize(&argc, &argv);
    int n = strtol(argv[1], nullptr, 10);
    int N = n * n;
    double h = 1.0 / n;

    double dx = 1.0, dy = 1.0;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;

    A.SetInterval(0, N);
    b.SetInterval(0, N);
    sol.SetInterval(0, N);

    for (int k = 0; k < N; ++k) {
        int i = k / n;
        int j = k % n;
        double x = h * i;
        double y = h * j;

        //boundary points
        if (i == 0 || j == 0 || i == n - 1 || j == n - 1) { 
            A[k][k] = 1;
            b[k] = g(x, y);
        } else {
            A[k][k] = 2.0 * (dx + dy);
            A[k][k + n] = -dx;
            A[k][k - n] = -dx;
            A[k][k - 1] = -dy;
            A[k][k + 1] = -dy;
            b[k] = f(x, y) * h * h;
        }
    }

    Solver S(Solver::INNER_ILU2);
    S.SetParameter("absolute_tolerance", "1e-12");
    S.SetParameter("verbosity", "2");
    S.SetParameter("relative_tolerance", "1e-7");
    S.SetMatrix(A);
    bool solved = S.Solve(b, sol);
    std::cout << "lin.it.: " << S.Iterations() << std::endl;
    if (!solved) {
        std::cout << "Linear solver failure!" << std::endl;
        std::cout << "Reason: " << S.ReturnReason() << std::endl;
    }

    std::vector<double> full_sol(N);
    for (int i = 0; i < N; ++i) {
        full_sol[i] = sol[i];
    }

    std::vector<double> exact_sol(N);
    for (int k = 0; k < N; ++k) {
        int i = k / n;
        int j = k % n;
        exact_sol[k] = u(i * h, j * h, dx, dy);
    }


    std::string res_fname = "./res.vtk";

    double l2_norm = 0., c_norm = 0;
    for (int i = 0; i < N; ++i) {
        l2_norm += (exact_sol[i] - sol[i]) * (exact_sol[i] - sol[i]);
        c_norm = std::max(c_norm, abs(exact_sol[i] - sol[i]));
    }
    l2_norm = sqrt(l2_norm) / N;
    std::ofstream norms("norms.txt");
    norms << "L2 norm: " << l2_norm << " C-norm: " << c_norm << std::endl;
    std::cout << "L2 norm: " << l2_norm << " C-norm: " << c_norm << std::endl;
    norms.close();

    Solver::Finalize();
}
