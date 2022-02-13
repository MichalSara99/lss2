#include "lss_double_sweep_solver_t.hpp"

#include "../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../boundaries/lss_neumann_boundary.hpp"
#include "../../boundaries/lss_robin_boundary.hpp"
#include "../../sparse_solvers/tridiagonal/double_sweep_solver/lss_double_sweep_solver.hpp"

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using lss_double_sweep_solver::double_sweep_solver;
using lss_utility::container_t;

void dss_bvp_dirichlet_bc_test()
{

    /*

    Solve BVP:

    u''(t) = - 2,

    where

    t \in (0, 1)
    u(0) = 0 ,  u(1) = 0


    Exact solution is

    u(t) = t(1-t)

    */
    std::cout << "=================================\n";
    std::cout << "Solving Boundary-value problem: \n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = t(1-t)\n";
    std::cout << "=================================\n";

    // discretization:
    std::size_t N{100};
    // step size:
    double h = 1.0 / static_cast<double>(N - 1);
    // upper,mid, and lower diagonal:
    container_t upper_diag(1.0, N);
    container_t diagonal(-2.0, N);
    container_t lower_diag(1.0, N);

    // right-hand side:
    container_t rhs(-2.0 * h * h, N);

    // boundary conditions:
    auto const &lower_ptr = std::make_shared<dirichlet_boundary_1d>([](double t) { return 0.0; });
    auto const &upper_ptr = std::make_shared<dirichlet_boundary_1d>([](double t) { return 0.0; });

    auto dss = std::make_shared<double_sweep_solver>(N);
    dss->set_diagonals(std::move(lower_diag), std::move(diagonal), std::move(upper_diag));
    dss->set_rhs(rhs);
    // get the solution:
    container_t solution(N);
    dss->solve(std::make_pair(lower_ptr, upper_ptr), solution);

    // exact value:
    auto exact = [](double x) { return x * (1.0 - x); };

    std::cout << "tp : FDM | Exact | Diff\n";
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j << ": " << solution[j] << " |  " << exact(j * h) << " | " << (solution[j] - exact(j * h))
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact(t * h)) < 1e-14, "element " << t << " must be the same.");
    }
}

void test_dss_bvp_dirichlet_bc()
{
    std::cout << "==================================================\n";
    std::cout << "=========== Double Sweep (Dirichlet BC) ==========\n";
    std::cout << "==================================================\n";

    dss_bvp_dirichlet_bc_test();

    std::cout << "==================================================\n";
}

void dss_bvp_robin_bc_test()
{

    /*

    Solve BVP:

    u''(t) = - 2,

    where

    t \in (0, 1)
    u(0) = 1 ,  u'(1) + u(1) = 0


    Exact solution is

    u(t) = -t*t + t + 1

    */
    std::cout << "=================================\n";
    std::cout << "Solving Boundary-value problem: \n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(1) + u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 1\n";
    std::cout << "=================================\n";

    // discretization:
    std::size_t N{100};
    // step size:
    auto const h = 1.0 / static_cast<double>(N - 1);
    // upper,mid, and lower diagonal:
    container_t upper_diag(1.0, N);
    container_t diagonal(-2.0, N);
    container_t lower_diag(1.0, N);

    // right-hand side:
    container_t rhs(-2.0 * h * h, N);

    // boundary conditions:
    auto const &lower_ptr = std::make_shared<dirichlet_boundary_1d>([](double t) { return 1.0; });
    auto const &upper_ptr =
        std::make_shared<robin_boundary_1d>([](double t) { return 1.0; }, [](double t) { return 0.0; });

    auto dss = std::make_shared<double_sweep_solver>(N);
    dss->set_diagonals(std::move(lower_diag), std::move(diagonal), std::move(upper_diag));
    dss->set_rhs(rhs);
    // get the solution:
    container_t solution(N);
    dss->solve(std::make_pair(lower_ptr, upper_ptr), solution);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + 1.0); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j << ": " << solution[j] << " |  " << exact(j * h) << " | " << (solution[j] - exact(j * h))
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact(t * h)) < 1e-13, "element " << t << " must be the same.");
    }
}

void test_dss_bvp_robin_bc()
{
    std::cout << "==================================================\n";
    std::cout << "=========== Double Sweep (Robin BC) ==============\n";
    std::cout << "==================================================\n";

    dss_bvp_robin_bc_test();

    std::cout << "==================================================\n";
}

void dss_bvp_mix_bc_test()
{

    /*

    Solve BVP:

    u''(t) = - 2,

    where

    t \in (0, 1)
    u'(0) - 1  = 0,  u'(1) + 2*u(1) = 0


    Exact solution is

    u(t) = -t*t + t + 0.5

    */
    std::cout << "=================================\n";
    std::cout << "Solving Boundary-value problem: \n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " u''(t) = -2, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,1>,\n";
    std::cout << " u'(0) - 1 = 0 \n";
    std::cout << " u'(1) + 2*u(1) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = -t*t + t + 0.5\n";
    std::cout << "=================================\n";

    // discretization:
    std::size_t N{1000};
    // step size:
    const auto h = 1.0 / static_cast<double>(N - 1);
    // upper,mid, and lower diagonal:
    container_t upper_diag(1.0, N);
    container_t diagonal(-2.0, N);
    container_t lower_diag(1.0, N);

    // right-hand side:
    container_t rhs(-2.0 * h * h, N);

    // boundary conditions:
    auto const &lower_ptr = std::make_shared<neumann_boundary_1d>([](double t) { return -1.0; });
    auto const &upper_ptr =
        std::make_shared<robin_boundary_1d>([](double t) { return 2.0; }, [](double t) { return 0.0; });
    // constriuct space range:
    auto dss = std::make_shared<double_sweep_solver>(N);
    dss->set_diagonals(std::move(lower_diag), std::move(diagonal), std::move(upper_diag));
    dss->set_rhs(rhs);
    // get the solution:
    container_t solution(N);
    dss->solve(std::make_pair(lower_ptr, upper_ptr), solution);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + 0.5); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j << ": " << solution[j] << " |  " << exact(j * h) << " | " << (solution[j] - exact(j * h))
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        LSS_ASSERT(std::abs(solution[t] - exact(t * h)) < 1e-13, "element " << t << " must be the same.");
    }
}

void test_dss_bvp_mix_bc()
{
    std::cout << "==================================================\n";
    std::cout << "=========== Double Sweep (Mix BC) ==============\n";
    std::cout << "==================================================\n";

    dss_bvp_mix_bc_test();

    std::cout << "==================================================\n";
}
