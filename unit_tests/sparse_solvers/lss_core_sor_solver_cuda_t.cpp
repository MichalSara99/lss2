#include "lss_core_sor_solver_cuda_t.hpp"

#include "../../sparse_solvers/core/cuda_sor_solver/lss_core_sor_solver_cuda.hpp"
#include <iostream>
#include <string>

using lss_core_sor_solver::core_sor_solver_cuda;
using lss_core_sor_solver::flat_matrix;
using lss_utility::container_t;

void device_sor_solver_basic_test()
{

    std::cout << "===================================\n";
    std::cout << "Solving sparse system of equations:\n";
    std::cout << "Ax = b,\n";
    std::cout << "where\n\n";
    std::cout << "     [10, -1, 2, 0] \n";
    std::cout << " A = [-1, 11, -1, 3] \n";
    std::cout << "     [2, -1, 10, -1] \n";
    std::cout << "	 [0, 3, -1, 8]\n";
    std::cout << "\n";
    std::cout << "b = [6, 25, -11, 15]^T \n";
    std::cout << "\n";

    /*

    Solving for x in

    Ax = b,
    where


    A =
    [10, -1, 2, 0]
    [-1, 11, -1, 3]
    [2, -1, 10, -1]
    [0, 3, -1, 8]

    b = [6, 25, -11, 15]^T

    */
    // size of the system:
    int const m = 4;
    // first create and populate the sparse matrix:
    flat_matrix fm(m, m);

    // populate the matrix:
    fm.emplace_back(0, 0, 10.0);
    fm.emplace_back(0, 1, -1.);
    fm.emplace_back(0, 2, 2.);

    fm.emplace_back(1, 0, -1.);
    fm.emplace_back(1, 1, 11.);
    fm.emplace_back(1, 2, -1.);
    fm.emplace_back(1, 3, 3.);

    fm.emplace_back(2, 0, 2.0);
    fm.emplace_back(2, 1, -1.0);
    fm.emplace_back(2, 2, 10.0);
    fm.emplace_back(2, 3, -1.);

    fm.emplace_back(3, 1, 3.0);
    fm.emplace_back(3, 2, -1.0);
    fm.emplace_back(3, 3, 8.);

    // lets use std::vector to populate vector b:
    container_t b = {6., 25.0, -11., 15.};

    // create sparse solver on DEVICE:
    core_sor_solver_cuda sor(m);

    // insert sparse matrix A and vector b:
    sor.set_flat_sparse_matrix(std::move(fm));
    sor.set_rhs(b);
    sor.set_omega(0.5);

    auto solution = sor.solve();
    std::cout << "Solution is: \n[";
    for (auto const &e : solution)
    {
        std::cout << e << " ";
    }
    std::cout << "]\n";
    container_t true_sol = {1.0, 2.0, -1.0, 1.0};
    for (std::size_t t = 0; t < true_sol.size(); ++t)
    {
        LSS_ASSERT(std::abs(true_sol[t] - solution[t]) < 0.000001, "element " << t << " must be the same");
    }
}

void test_device_sor_solver_basic()
{
    std::cout << "============================================================\n";
    std::cout << "=============== Initialise Flat-Raw-Matrix =================\n";
    std::cout << "============================================================\n";

    device_sor_solver_basic_test();

    std::cout << "============================================================\n";
}

void device_sor_solver_ode_test()
{

    std::cout << "=================================\n";
    std::cout << " Using SOR algorithm to \n";
    std::cout << " solve Boundary Value Problem: \n\n";
    std::cout << " type: " << typeid(double).name() << "\n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = t(1-t)\n";
    std::cout << "=================================\n";

    // discretization:
    // t_0,t_1,t_2,...,t_20
    int const N = 20;
    // step size:
    double h = 1.0 / static_cast<double>(N);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:

    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fm(m, m);

    // populate the matrix:
    fm.emplace_back(0, 0, -2.0);
    fm.emplace_back(0, 1, 1.0);
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fm.emplace_back(t, t - 1, 1.0);
        fm.emplace_back(t, t, -2.0);
        fm.emplace_back(t, t + 1, 1.0);
    }
    fm.emplace_back(m - 1, m - 2, 1.0);
    fm.emplace_back(m - 1, m - 1, -2.0);

    // lets use std::vector to populate vector b:
    container_t b(-2.0 * h * h, m);
    // set the Dirichlet boundary conditions:
    auto const left = 0.0;
    auto const right = 0.0;
    b[0] = b[0] - left;
    b[b.size() - 1] = b[b.size() - 1] - right;

    // create sor  solver:
    core_sor_solver_cuda sor(m);

    // insert sparse matrix A and vector b:
    sor.set_flat_sparse_matrix(std::move(fm));
    sor.set_rhs(b);
    sor.set_omega(0.5);

    container_t solution(m);
    sor.solve(solution);

    // exact value:
    auto exact = [](double x) { return x * (1.0 - x); };

    std::cout << "tp : FDM | Exact\n";
    std::cout << "t_" << 0 << ": " << left << " |  " << exact(0) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << '\n';
    }
    std::cout << "t_" << N << ": " << right << " |  " << exact(N * h) << '\n';
    // TEST:
    LSS_ASSERT(std::abs(left - exact(0)) < 0.000001, "most left-most element  must be the same");
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        LSS_ASSERT(std::abs(exact((j + 1) * h) - solution[j]) < 0.000001, "element " << j << " must be the same");
    }
    LSS_ASSERT(std::abs(right - exact(N * h)) < 0.000001, "most right-most element  must be the same");
}

void test_device_sor_solver_ode()
{
    std::cout << "============================================================\n";
    std::cout << "=========================== BVP ============================\n";
    std::cout << "============================================================\n";

    device_sor_solver_ode_test();

    std::cout << "============================================================\n";
}
