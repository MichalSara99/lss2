#include "lss_karawia_solver_t.hpp"

#include "../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../sparse_solvers/pentadiagonal/karawia_solver/lss_karawia_solver.hpp"

using lss_boundary::dirichlet_boundary_1d;
using lss_karawia_solver::karawia_solver;
using lss_utility::container_t;

void karawia_bvp_dirichlet_bc_test()
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
    container_t uppest_diag(-1.0, N);
    container_t upper_diag(16.0, N);
    container_t diagonal(-30.0, N);
    container_t lower_diag(16.0, N);
    container_t lowest_diag(-1.0, N);

    // right-hand side:
    container_t rhs(-2.0 * 12.0 * h * h, N);

    // exact value:
    auto exact = [](double x) { return x * (1.0 - x); };

    // boundary conditions:
    auto const &lower_ptr = std::make_shared<dirichlet_boundary_1d>([](double t) { return 0.0; });
    auto const &upper_ptr = std::make_shared<dirichlet_boundary_1d>([](double t) { return 0.0; });
    auto const &other_lower_ptr = std::make_shared<dirichlet_boundary_1d>([&](double t) { return exact(h); });
    auto const &other_upper_ptr = std::make_shared<dirichlet_boundary_1d>([&](double t) { return exact(1.0 - h); });

    auto dss = std::make_shared<karawia_solver>(N);
    dss->set_diagonals(std::move(lowest_diag), std::move(lower_diag), std::move(diagonal), std::move(upper_diag),
                       std::move(uppest_diag));
    dss->set_rhs(rhs);
    // get the solution:
    container_t solution(N);
    //
    auto const &boundary = std::make_pair(lower_ptr, upper_ptr);
    auto const &other_boundary = std::make_pair(other_lower_ptr, other_upper_ptr);
    dss->solve(boundary, other_boundary, solution);

    std::cout << "tp : FDM | Exact | Diff\n";
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        std::cout << "t_" << j << ": " << solution[j] << " |  " << exact(j * h) << " | " << (solution[j] - exact(j * h))
                  << '\n';
    }
}

void test_karawia_bvp_dirichlet_bc()
{
    std::cout << "==================================================\n";
    std::cout << "============== Karawia (Dirichlet BC) ============\n";
    std::cout << "==================================================\n";

    karawia_bvp_dirichlet_bc_test();

    std::cout << "==================================================\n";
}
