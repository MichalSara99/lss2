#include "lss_core_cuda_solver_t.hpp"

#include "../../sparse_solvers/core/cuda_solver/lss_core_cuda_solver.hpp"
#include "../../sparse_solvers/core/cuda_solver/lss_core_cuda_solver_policy.hpp"
#include <iostream>
#include <string>

using lss_core_cuda_solver::factorization_enum;
using lss_core_cuda_solver::flat_matrix;
using lss_core_cuda_solver::memory_space_enum;
using lss_core_cuda_solver::real_sparse_solver_cuda;
using lss_core_cuda_solver_policy::sparse_solver_device_qr;
using lss_core_cuda_solver_policy::sparse_solver_host_qr;
using lss_utility::container_t;

void device_sparse_qr_test()
{
    std::cout << "===================================\n";
    std::cout << "Solving sparse system of equations:\n";
    std::cout << "Ax = b,\n";
    std::cout << "where\n\n";
    std::cout << "    [  1.0  2.0  0.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  3.0  4.0  5.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  0.0  6.0  7.0  8.0  0.0  0.0 ] \n";
    std::cout << "A = [  0.0  0.0  9.0 10.0 11.0  0.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0 12.0 13.0 14.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0  0.0 15.0 16.0 ] \n";
    std::cout << "\n";
    std::cout << "b = [  0.0  2.0  4.0  6.0  8.0  10.0 ]^T \n";
    std::cout << "\n";

    /*

    Solving for x in

    Ax = b,
    where


    A =
    [  1.0  2.0  0.0  0.0  0.0  0.0 ]
    [  3.0  4.0  5.0  0.0  0.0  0.0 ]
    [  0.0  6.0  7.0  8.0  0.0  0.0 ]
    [  0.0  0.0  9.0 10.0 11.0  0.0 ]
    [  0.0  0.0  0.0 12.0 13.0 14.0 ]
    [  0.0  0.0  0.0  0.0 15.0 16.0 ]

    b = [0.0,2.0,4.0,6.0,8.0,10.0]^T

    */

    // size of the system:
    int const m = 6;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (1.0));
    fsm.emplace_back(0, 1, (2.0));
    fsm.emplace_back(1, 0, (3.0));
    fsm.emplace_back(1, 1, (4.0));
    fsm.emplace_back(1, 2, (5.0));
    fsm.emplace_back(2, 1, (6.0));
    fsm.emplace_back(2, 2, (7.0));
    fsm.emplace_back(2, 3, (8.0));
    fsm.emplace_back(3, 2, (9.0));
    fsm.emplace_back(3, 3, (10.0));
    fsm.emplace_back(3, 4, (11.0));
    fsm.emplace_back(4, 3, (12.0));
    fsm.emplace_back(4, 4, (13.0));
    fsm.emplace_back(4, 5, (14.0));
    fsm.emplace_back(5, 4, (15.0));
    fsm.emplace_back(5, 5, (16.0));

    // lets use std::vector to populate vector b:
    container_t b = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

    // true solution is:
    container_t true_sol = {-0.382378, 0.191189, 0.476476, -0.060308, 0.210436, 0.427716};

    // create sparse solver on DEVICE:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    auto solution = rss.solve();
    std::cout << "[";
    for (std::size_t t = 0; t < m; ++t)
    {
        std::cout << solution[t] << " ";
    }
    std::cout << "]\n";
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - true_sol[t]) < 0.00001, "element " << t << " must equal");
    }
}

extern void device_sparse_qr_pointer_test()
{

    std::cout << "===================================\n";
    std::cout << "Solving sparse system of equations:\n";
    std::cout << "Ax = b,\n";
    std::cout << "where\n\n";
    std::cout << "    [  1.0  2.0  0.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  3.0  4.0  5.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  0.0  6.0  7.0  8.0  0.0  0.0 ] \n";
    std::cout << "A = [  0.0  0.0  9.0 10.0 11.0  0.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0 12.0 13.0 14.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0  0.0 15.0 16.0 ] \n";
    std::cout << "\n";
    std::cout << "b = [  0.0  2.0  4.0  6.0  8.0  10.0 ]^T \n";
    std::cout << "\n";

    /*

    Solving for x in

    Ax = b,
    where


    A =
    [  1.0  2.0  0.0  0.0  0.0  0.0 ]
    [  3.0  4.0  5.0  0.0  0.0  0.0 ]
    [  0.0  6.0  7.0  8.0  0.0  0.0 ]
    [  0.0  0.0  9.0 10.0 11.0  0.0 ]
    [  0.0  0.0  0.0 12.0 13.0 14.0 ]
    [  0.0  0.0  0.0  0.0 15.0 16.0 ]

    b = [0.0,2.0,4.0,6.0,8.0,10.0]^T

    */

    // size of the system:
    int const m = 6;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (1.0));
    fsm.emplace_back(0, 1, (2.0));
    fsm.emplace_back(1, 0, (3.0));
    fsm.emplace_back(1, 1, (4.0));
    fsm.emplace_back(1, 2, (5.0));
    fsm.emplace_back(2, 1, (6.0));
    fsm.emplace_back(2, 2, (7.0));
    fsm.emplace_back(2, 3, (8.0));
    fsm.emplace_back(3, 2, (9.0));
    fsm.emplace_back(3, 3, (10.0));
    fsm.emplace_back(3, 4, (11.0));
    fsm.emplace_back(4, 3, (12.0));
    fsm.emplace_back(4, 4, (13.0));
    fsm.emplace_back(4, 5, (14.0));
    fsm.emplace_back(5, 4, (15.0));
    fsm.emplace_back(5, 5, (16.0));

    // lets use std::vector to populate vector b:
    container_t b = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

    // true solution is:
    container_t true_sol = {-0.382378, 0.191189, 0.476476, -0.060308, 0.210436, 0.427716};

    // create sparse solver on DEVICE:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution);

    std::cout << "Solution is: \n[";
    for (std::size_t t = 0; t < m; ++t)
    {
        std::cout << solution[t] << " ";
    }
    std::cout << "]\n";
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - true_sol[t]) < 0.00001, "element " << t << " must equal");
    }
}

extern void test_device_sparse_qr()
{
    std::cout << "==================================================\n";
    std::cout << "=========== Sparse QR Solver - DEVICE ============\n";
    std::cout << "==================================================\n";

    device_sparse_qr_test();
    device_sparse_qr_pointer_test();

    std::cout << "==================================================\n";
}

extern void host_sparse_qr_test()
{

    std::cout << "===================================\n";
    std::cout << "Solving sparse system of equations:\n";
    std::cout << "Ax = b,\n";
    std::cout << "where\n\n";
    std::cout << "    [  1.0  2.0  0.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  3.0  4.0  5.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  0.0  6.0  7.0  8.0  0.0  0.0 ] \n";
    std::cout << "A = [  0.0  0.0  9.0 10.0 11.0  0.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0 12.0 13.0 14.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0  0.0 15.0 16.0 ] \n";
    std::cout << "\n";
    std::cout << "b = [  0.0  2.0  4.0  6.0  8.0  10.0 ]^T \n";
    std::cout << "\n";

    /*

    Solving for x in

    Ax = b,
    where


    A =
    [  1.0  2.0  0.0  0.0  0.0  0.0 ]
    [  3.0  4.0  5.0  0.0  0.0  0.0 ]
    [  0.0  6.0  7.0  8.0  0.0  0.0 ]
    [  0.0  0.0  9.0 10.0 11.0  0.0 ]
    [  0.0  0.0  0.0 12.0 13.0 14.0 ]
    [  0.0  0.0  0.0  0.0 15.0 16.0 ]

    b = [0.0,2.0,4.0,6.0,8.0,10.0]^T

    */

    // size of the system:
    int const m = 6;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (1.0));
    fsm.emplace_back(0, 1, (2.0));
    fsm.emplace_back(1, 0, (3.0));
    fsm.emplace_back(1, 1, (4.0));
    fsm.emplace_back(1, 2, (5.0));
    fsm.emplace_back(2, 1, (6.0));
    fsm.emplace_back(2, 2, (7.0));
    fsm.emplace_back(2, 3, (8.0));
    fsm.emplace_back(3, 2, (9.0));
    fsm.emplace_back(3, 3, (10.0));
    fsm.emplace_back(3, 4, (11.0));
    fsm.emplace_back(4, 3, (12.0));
    fsm.emplace_back(4, 4, (13.0));
    fsm.emplace_back(4, 5, (14.0));
    fsm.emplace_back(5, 4, (15.0));
    fsm.emplace_back(5, 5, (16.0));

    // lets use std::vector to populate vector b:
    container_t b = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

    // true solution is:
    container_t true_sol = {-0.382378, 0.191189, 0.476476, -0.060308, 0.210436, 0.427716};

    // create sparse solver on DEVICE:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    auto solution = rss.solve();
    std::cout << "Solution is: \n[";
    for (auto const &e : solution)
    {
        std::cout << e << " ";
    }
    std::cout << "]\n";
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - true_sol[t]) < 0.00001, "element " << t << " must equal");
    }
}

extern void host_sparse_qr_pointer_test()
{

    std::cout << "===================================\n";
    std::cout << "Solving sparse system of equations:\n";
    std::cout << "Ax = b,\n";
    std::cout << "where\n\n";
    std::cout << "    [  1.0  2.0  0.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  3.0  4.0  5.0  0.0  0.0  0.0 ] \n";
    std::cout << "    [  0.0  6.0  7.0  8.0  0.0  0.0 ] \n";
    std::cout << "A = [  0.0  0.0  9.0 10.0 11.0  0.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0 12.0 13.0 14.0 ] \n";
    std::cout << "    [  0.0  0.0  0.0  0.0 15.0 16.0 ] \n";
    std::cout << "\n";
    std::cout << "b = [  0.0  2.0  4.0  6.0  8.0  10.0 ]^T \n";
    std::cout << "\n";

    /*

    Solving for x in

    Ax = b,
    where


    A =
    [  1.0  2.0  0.0  0.0  0.0  0.0 ]
    [  3.0  4.0  5.0  0.0  0.0  0.0 ]
    [  0.0  6.0  7.0  8.0  0.0  0.0 ]
    [  0.0  0.0  9.0 10.0 11.0  0.0 ]
    [  0.0  0.0  0.0 12.0 13.0 14.0 ]
    [  0.0  0.0  0.0  0.0 15.0 16.0 ]

    b = [0.0,2.0,4.0,6.0,8.0,10.0]^T

    */

    // size of the system:
    int const m = 6;

    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (1.0));
    fsm.emplace_back(0, 1, (2.0));
    fsm.emplace_back(1, 0, (3.0));
    fsm.emplace_back(1, 1, (4.0));
    fsm.emplace_back(1, 2, (5.0));
    fsm.emplace_back(2, 1, (6.0));
    fsm.emplace_back(2, 2, (7.0));
    fsm.emplace_back(2, 3, (8.0));
    fsm.emplace_back(3, 2, (9.0));
    fsm.emplace_back(3, 3, (10.0));
    fsm.emplace_back(3, 4, (11.0));
    fsm.emplace_back(4, 3, (12.0));
    fsm.emplace_back(4, 4, (13.0));
    fsm.emplace_back(4, 5, (14.0));
    fsm.emplace_back(5, 4, (15.0));
    fsm.emplace_back(5, 5, (16.0));

    // lets use std::vector to populate vector b:
    container_t b = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};

    // true solution is:
    container_t true_sol = {-0.382378, 0.191189, 0.476476, -0.060308, 0.210436, 0.427716};

    // create sparse solver on DEVICE:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    auto solution = rss.solve();
    std::cout << "Solution is: \n[";
    for (auto const &e : solution)
    {
        std::cout << e << " ";
    }
    std::cout << "]\n";
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - true_sol[t]) < 0.00001, "element " << t << " must equal");
    }
}

extern void test_host_sparse_qr_test()
{
    std::cout << "==================================================\n";
    std::cout << "=========== Sparse QR Solver - HOST ==============\n";
    std::cout << "==================================================\n";

    host_sparse_qr_test();
    host_sparse_qr_pointer_test();

    std::cout << "==================================================\n";
}

extern void host_bvp_dirichlet_bc_qr_test()
{

    std::cout << "=================================\n";
    std::cout << " Using QR decomposition to \n";
    std::cout << " solve Boundary Value Problem: \n\n";
    std::cout << " type: double					\n\n";
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
    const double h = (1.0) / (N);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (-2.0));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, (-2.0));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // set the Dirichlet boundary conditions:
    double left = (0.0);
    double right = (0.0);
    b[0] = b[0] - left;
    b[b.size() - 1] = b[b.size() - 1] - right;

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution);

    // exact value:
    auto exact = [](double x) { return x * ((1.0) - x); };

    std::cout << "tp : FDM | Exact\n";
    std::cout << "t_" << 0 << ": " << left << " |  " << exact(0) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << '\n';
    }
    std::cout << "t_" << N << ": " << right << " |  " << exact(N * h) << '\n';

    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact((t + 1) * h)) < 0.00001, "element " << t << " must equal");
    }
}

extern void host_bvp_dirichlet_bc_lu_test()
{

    std::cout << "=================================\n";
    std::cout << " Using LU decomposition to \n";
    std::cout << " solve Boundary Value Problem: \n\n";
    std::cout << " type: double					\n\n";
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
    const double h = (1.0) / (N);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (-2.0));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, (-2.0));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // set the Dirichlet boundary conditions:
    const double left = (0.0);
    const double right = (0.0);
    b[0] = b[0] - left;
    b[b.size() - 1] = b[b.size() - 1] - right;

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution, factorization_enum::LUMethod);

    // exact value:
    auto exact = [](double x) { return x * ((1.0) - x); };

    std::cout << "tp : FDM | Exact\n";
    std::cout << "t_" << 0 << ": " << left << " |  " << exact(0) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << '\n';
    }
    std::cout << "t_" << N << ": " << right << " |  " << exact(N * h) << '\n';
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact((t + 1) * h)) < 0.00001, "element " << t << " must equal");
    }
}

extern void host_bvp_dirichlet_bc_cholesky_test()
{

    std::cout << "=================================\n";
    std::cout << " Using Cholesky decomposition to \n";
    std::cout << " solve Boundary Value Problem: \n\n";
    std::cout << " type: double					\n\n";
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
    double h = (1.0) / (N);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (-2.0));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, (-2.0));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // set the Dirichlet boundary conditions:
    double left = (0.0);
    double right = (0.0);
    b[0] = b[0] - left;
    b[b.size() - 1] = b[b.size() - 1] - right;

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution, factorization_enum::CholeskyMethod);

    // exact value:
    auto exact = [](double x) { return x * ((1.0) - x); };

    std::cout << "tp : FDM | Exact\n";
    std::cout << "t_" << 0 << ": " << left << " |  " << exact(0) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << '\n';
    }
    std::cout << "t_" << N << ": " << right << " |  " << exact(N * h) << '\n';
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact((t + 1) * h)) < 0.00001, "element " << t << " must equal");
    }
}

extern void test_dirichlet_bc_bvp_on_host()
{
    std::cout << "==================================================\n";
    std::cout << "============ Dirichlet BC BVP - HOST ´´============\n";
    std::cout << "==================================================\n";

    host_bvp_dirichlet_bc_qr_test();
    host_bvp_dirichlet_bc_lu_test();
    host_bvp_dirichlet_bc_cholesky_test();

    std::cout << "==================================================\n";
}

extern void device_bvp_dirichlet_bc_qr_test()
{

    std::cout << "=================================\n";
    std::cout << " Using QR decomposition to \n";
    std::cout << " solve Boundary Value Problem: \n\n";
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
    double h = (1.0) / (N);
    // set number of columns and rows:
    // because we already know the boundary values
    // at timepoints t_0 and t_20:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (-2.0));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, (-2.0));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // set the Dirichlet boundary conditions:
    double left = (0.0);
    double right = (0.0);
    b[0] = b[0] - left;
    b[b.size() - 1] = b[b.size() - 1] - right;

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution);

    // exact value:
    auto exact = [](double x) { return x * ((1.0) - x); };

    std::cout << "tp : FDM | Exact\n";
    std::cout << "t_" << 0 << ": " << left << " |  " << exact(0) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << '\n';
    }
    std::cout << "t_" << N << ": " << right << " |  " << exact(N * h) << '\n';
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact((t + 1) * h)) < 0.00001, "element " << t << " must equal");
    }
}

extern void device_bvp_dirichlet_bc_cholesky_test()
{

    std::cout << "=================================\n";
    std::cout << " Using Cholesky decomposition to \n";
    std::cout << " solve Boundary Value Problem: \n\n";
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
    double h = (1.0) / (N);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, (-2.0));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, (-2.0));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // set the Dirichlet boundary conditions:
    double left = (0.0);
    double right = (0.0);
    b[0] = b[0] - left;
    b[b.size() - 1] = b[b.size() - 1] - right;

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution, factorization_enum::CholeskyMethod);

    // exact value:
    auto exact = [](double x) { return x * ((1.0) - x); };

    std::cout << "tp : FDM | Exact\n";
    std::cout << "t_" << 0 << ": " << left << " |  " << exact(0) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << '\n';
    }
    std::cout << "t_" << N << ": " << right << " |  " << exact(N * h) << '\n';
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact((t + 1) * h)) < 0.00001, "element " << t << " must equal");
    }
}

extern void test_dirichlet_bc_bvp_on_device()
{
    std::cout << "==================================================\n";
    std::cout << "============ Dirichlet BC BVP - DEVICE ===========\n";
    std::cout << "==================================================\n";

    device_bvp_dirichlet_bc_qr_test();
    device_bvp_dirichlet_bc_cholesky_test();

    std::cout << "==================================================\n";
}

extern void host_bvp_robin_bc_qr_test()
{

    std::cout << "=================================\n";
    std::cout << " Using QR decomposition to \n";
    std::cout << " solve Boundary-value problem: \n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(1) + u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 1\n";
    std::cout << "=================================\n";

    // discretization:
    // t_0,t_1,t_2,...,t_20
    int const N = 100;
    // step size:
    double h = (1.0) / static_cast<double>(N);
    // set the Robin boundary conditions:
    double alpha = (.00);
    double phi = (1.0);
    double beta = ((2.0) + h) / ((2.0) - h);
    double psi = (0.0);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, ((alpha * 1.0 - 2.0)));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, ((-2.0 + (1.0 / beta))));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // adjust first and last elements due to the Robin BC
    b[0] = b[0] - (1.0) * phi;
    b[b.size() - 1] = b[b.size() - 1] + psi * ((1.0) / beta);

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + (1.0)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    std::cout << "t_" << 0 << ": " << (alpha * solution[0] + phi) << " |  " << exact(0) << " | "
              << ((alpha * solution[0] + phi) - exact(0)) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << " | "
                  << (solution[j] - exact((j + 1) * h)) << '\n';
    }
    std::cout << "t_" << N << ": " << ((solution[m - 1] - psi) / beta) << " |  " << exact(N * h) << " | "
              << (((solution[m - 1] - psi) / beta) - exact(N * h)) << '\n';
}

extern void host_bvp_robin_bc_lu_test()
{

    std::cout << "=================================\n";
    std::cout << " Using LU decomposition to \n";
    std::cout << " solve Boundary-value problem: \n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(1) + u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 1\n";
    std::cout << "=================================\n";

    // discretization:
    // t_0,t_1,t_2,...,t_20
    int const N = 100;
    // step size:
    double h = (1.0) / (N);
    // set the Robin boundary conditions:
    double alpha = (.00);
    double phi = (1.0);
    double beta = ((2.0) + h) / ((2.0) - h);
    double psi = (0.0);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, ((alpha * 1.0 - 2.0)));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, ((-2.0 + (1.0 / beta))));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // adjust first and last elements due to the Robin BC
    b[0] = b[0] - (1.0) * phi;
    b[b.size() - 1] = b[b.size() - 1] + psi * ((1.0) / beta);

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution, factorization_enum::LUMethod);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + (1.0)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    std::cout << "t_" << 0 << ": " << (alpha * solution[0] + phi) << " |  " << exact(0) << " | "
              << ((alpha * solution[0] + phi) - exact(0)) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << " | "
                  << (solution[j] - exact((j + 1) * h)) << '\n';
    }
    std::cout << "t_" << N << ": " << ((solution[m - 1] - psi) / beta) << " |  " << exact(N * h) << " | "
              << (((solution[m - 1] - psi) / beta) - exact(N * h)) << '\n';
}

extern void host_bvp_robin_bc_cholesky_test()
{

    std::cout << "=================================\n";
    std::cout << " Using Cholesky decomposition to \n";
    std::cout << " solve Boundary-value problem: \n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(1) + u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 1\n";
    std::cout << "=================================\n";

    // discretization:
    // t_0,t_1,t_2,...,t_20
    int const N = 100;
    // step size:
    double h = (1.0) / (N);
    // set the Robin boundary conditions:
    double alpha = (.00);
    double phi = (1.0);
    double beta = ((2.0) + h) / ((2.0) - h);
    double psi = (0.0);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, ((alpha * 1.0 - 2.0)));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, ((-2.0 + (1.0 / beta))));

    // lets use std::vector to populate vector b:
    container_t b((-2.0) * h * h, m);
    // adjust first and last elements due to the Robin BC
    b[0] = b[0] - (1.0) * phi;
    b[b.size() - 1] = b[b.size() - 1] + psi * ((1.0) / beta);

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Host> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution, factorization_enum::CholeskyMethod);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + (1.0)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    std::cout << "t_" << 0 << ": " << (alpha * solution[0] + phi) << " |  " << exact(0) << " | "
              << ((alpha * solution[0] + phi) - exact(0)) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << " | "
                  << (solution[j] - exact((j + 1) * h)) << '\n';
    }
    std::cout << "t_" << N << ": " << ((solution[m - 1] - psi) / beta) << " |  " << exact(N * h) << " | "
              << (((solution[m - 1] - psi) / beta) - exact(N * h)) << '\n';
}

extern void test_robin_bc_bvp_on_host()
{
    std::cout << "==================================================\n";
    std::cout << "============ Robin BC BVP - HOST =================\n";
    std::cout << "==================================================\n";

    host_bvp_robin_bc_qr_test();
    host_bvp_robin_bc_lu_test();
    // host_bvp_robin_bc_cholesky_test();

    std::cout << "==================================================\n";
}

extern void device_bvp_robin_bc_qr_test()
{

    std::cout << "=================================\n";
    std::cout << " Using QR decomposition to \n";
    std::cout << " solve Boundary-value problem: \n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(1) + u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 1\n";
    std::cout << "=================================\n";

    // discretization:
    // t_0,t_1,t_2,...,t_20
    int const N = 100;
    // step size:
    double h = (1.0) / static_cast<double>(N);
    // set the Robin boundary conditions:
    double alpha = (.00);
    double phi = (1.0);
    double beta = ((2.0) + h) / ((2.0) - h);
    double psi = (0.0);
    // set number of columns and rows:
    // because we already know the boundary values
    // at timepoints t_0 and t_20:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, ((alpha * 1.0 - 2.0)));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, ((-2.0 + (1.0 / beta))));

    // lets use std::vector to populate vector b:
    container_t b(static_cast<double>(-2.0) * h * h, m);
    // adjust first and last elements due to the Robin BC
    b[0] = b[0] - static_cast<double>(1.0) * phi;
    b[b.size() - 1] = b[b.size() - 1] + psi * (static_cast<double>(1.0) / beta);

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + static_cast<double>(1.0)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    std::cout << "t_" << 0 << ": " << (alpha * solution[0] + phi) << " |  " << exact(0) << " | "
              << ((alpha * solution[0] + phi) - exact(0)) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << " | "
                  << (solution[j] - exact((j + 1) * h)) << '\n';
    }
    std::cout << "t_" << N << ": " << ((solution[m - 1] - psi) / beta) << " |  " << exact(N * h) << " | "
              << (((solution[m - 1] - psi) / beta) - exact(N * h)) << '\n';
}

extern void device_bvp_robin_bc_cholesky_test()
{

    std::cout << "=================================\n";
    std::cout << " Using Cholesky decomposition to \n";
    std::cout << " solve Boundary-value problem: \n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(1) + u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 1\n";
    std::cout << "=================================\n";

    // discretization:
    // t_0,t_1,t_2,...,t_20
    int const N = 100;
    // step size:
    double h = static_cast<double>(1.0) / static_cast<double>(N);
    // set the Robin boundary conditions:
    double alpha = (.00);
    double phi = (1.0);
    double beta = ((2.0) + h) / ((2.0) - h);
    double psi = (0.0);
    // set number of columns and rows:
    // because we already know the boundary values
    // at t_0 = 0 and t_20 = 0:
    int const m = N - 1;
    // first create and populate the sparse matrix:
    flat_matrix fsm(m, m);
    // populate the matrix:
    fsm.emplace_back(0, 0, ((alpha * 1.0 - 2.0)));
    fsm.emplace_back(0, 1, (1.0));
    for (std::size_t t = 1; t < m - 1; ++t)
    {
        fsm.emplace_back(t, t - 1, (1.0));
        fsm.emplace_back(t, t, (-2.0));
        fsm.emplace_back(t, t + 1, (1.0));
    }
    fsm.emplace_back(m - 1, m - 2, (1.0));
    fsm.emplace_back(m - 1, m - 1, ((-2.0 + (1.0 / beta))));

    // lets use std::vector to populate vector b:
    container_t b(static_cast<double>(-2.0) * h * h, m);
    // adjust first and last elements due to the Robin BC
    b[0] = b[0] - static_cast<double>(1.0) * phi;
    b[b.size() - 1] = b[b.size() - 1] + psi * (static_cast<double>(1.0) / beta);

    // create sparse solver on HOST:
    real_sparse_solver_cuda<memory_space_enum::Device> rss(m);

    // because we used default cstor we need to call initialize
    rss.initialize(m);

    // insert sparse matrix A and vector b:
    rss.set_flat_sparse_matrix(std::move(fsm));
    rss.set_rhs(b);

    container_t solution(m);
    rss.solve(solution, factorization_enum::CholeskyMethod);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + static_cast<double>(1.0)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    std::cout << "t_" << 0 << ": " << (alpha * solution[0] + phi) << " |  " << exact(0) << " | "
              << ((alpha * solution[0] + phi) - exact(0)) << '\n';
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j + 1 << ": " << solution[j] << " |  " << exact((j + 1) * h) << " | "
                  << (solution[j] - exact((j + 1) * h)) << '\n';
    }
    std::cout << "t_" << N << ": " << ((solution[m - 1] - psi) / beta) << " |  " << exact(N * h) << " | "
              << (((solution[m - 1] - psi) / beta) - exact(N * h)) << '\n';
}

extern void test_robin_bc_bvp_on_device()
{
    std::cout << "==================================================\n";
    std::cout << "============ Robin BC BVP - DEVICE ===============\n";
    std::cout << "==================================================\n";

    device_bvp_robin_bc_qr_test();
    // device_bvp_robin_bc_cholesky_test();

    std::cout << "==================================================\n";
}
