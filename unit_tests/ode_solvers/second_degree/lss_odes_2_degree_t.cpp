#include "lss_odes_2_degree_t.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../ode_solvers/second_degree/lss_ode_equation.hpp"

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using lss_enumerations::grid_enum;
using lss_grids::grid_1d;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d;
using lss_grids::grid_transform_config_1d;
using lss_ode_solvers::ode_coefficient_data_config;
using lss_ode_solvers::ode_data_config;
using lss_ode_solvers::ode_discretization_config;
using lss_ode_solvers::ode_nonhom_data_config;
using lss_ode_solvers::default_ode_solver_configs::dev_cusolver_qr_solver_config_ptr;
using lss_ode_solvers::implicit_solvers::ode_equation;
using lss_utility::container_t;
using lss_utility::range;

void impl_simple_ode_dirichlet_bc_cuda_solver_device_qr()
{
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

    // number of space subdivisions:
    std::size_t Sd{100};
    // space range:
    auto const space_range = std::make_shared<range>(0.0, 1.0);
    // discretization config:
    auto const discretization_ptr = std::make_shared<ode_discretization_config>(space_range, Sd);
    // coeffs:
    auto a = [](double x) { return 0.0; };
    auto b = [](double x) { return 0.0; };
    auto const ode_coeffs_data_ptr = std::make_shared<ode_coefficient_data_config>(a, b);
    // nonhom data:
    auto two = [](double x) { return -2.0; };
    auto const ode_nonhom_data_ptr = std::make_shared<ode_nonhom_data_config>(two);
    // ode data config:
    auto const ode_data_ptr = std::make_shared<ode_data_config>(ode_coeffs_data_ptr, ode_nonhom_data_ptr);
    // boundary conditions:
    auto const &dirichlet = 0.0;
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize ode solver
    ode_equation odesolver(ode_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                           dev_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    odesolver.solve(solution);

    // exact value:
    auto exact = [](double x) { return x * (1.0 - x); };
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 7e-5, "element at " << t << "must be appr same.");
    }
}

void test_impl_simple_ode_dirichlet_bc_cuda_solver_device()
{
    std::cout << "==================================================\n";
    std::cout << "==================== (Dir-Dir BC) ================\n";
    std::cout << "==================================================\n";

    impl_simple_ode_dirichlet_bc_cuda_solver_device_qr();

    std::cout << "==================================================\n";
}

void impl_simple_ode_dirichlet_neumann_bc_cuda_solver_device_qr()
{

    std::cout << "=================================\n";
    std::cout << "Solving Boundary-value problem: \n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " u''(t) = 6*t, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,2>,\n";
    std::cout << " u(0) = 1 \n";
    std::cout << " u'(2) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = t*t*t - 12*t + 1\n";
    std::cout << "=================================\n";

    // number of space subdivisions:
    std::size_t Sd{100};
    // space range:
    auto const space_range = std::make_shared<range>(0.0, 2.0);
    // discretization config:
    auto const discretization_ptr = std::make_shared<ode_discretization_config>(space_range, Sd);
    // coeffs:
    auto a = [](double x) { return 0.0; };
    auto b = [](double x) { return 0.0; };
    auto const ode_coeffs_data_ptr = std::make_shared<ode_coefficient_data_config>(a, b);
    // nonhom data:
    auto two = [](double x) { return 6.0 * x; };
    auto const ode_nonhom_data_ptr = std::make_shared<ode_nonhom_data_config>(two);
    // ode data config:
    auto const ode_data_ptr = std::make_shared<ode_data_config>(ode_coeffs_data_ptr, ode_nonhom_data_ptr);
    // boundary conditions:
    auto const dirichlet = 1.0;
    auto const neumann = 0.0;
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize ode solver
    ode_equation odesolver(ode_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                           dev_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    odesolver.solve(solution);

    // exact value:
    auto exact = [](double x) { return (x * x * x - 12.0 * x + 1.0); };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-2, "element at " << t << "must be appr same.");
    }
}

void test_impl_simple_ode_dirichlet_neumann_bc_cuda_solver_device()
{
    std::cout << "==================================================\n";
    std::cout << "================= (Dir-Neu BC) ===================\n";
    std::cout << "==================================================\n";

    impl_simple_ode_dirichlet_neumann_bc_cuda_solver_device_qr();

    std::cout << "==================================================\n";
}

void impl_simple_ode_dirichlet_robin_bc_cuda_solver_device_qr()
{

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

    // number of space subdivisions:
    std::size_t Sd{100};
    // space range:
    auto const space_range = std::make_shared<range>(0.0, 1.0);
    // discretization config:
    auto const discretization_ptr = std::make_shared<ode_discretization_config>(space_range, Sd);
    // coeffs:
    auto a = [](double x) { return 0.0; };
    auto b = [](double x) { return 0.0; };
    auto const ode_coeffs_data_ptr = std::make_shared<ode_coefficient_data_config>(a, b);
    // nonhom data:
    auto two = [](double x) { return -2.0; };
    auto const ode_nonhom_data_ptr = std::make_shared<ode_nonhom_data_config>(two);
    // ode data config:
    auto const ode_data_ptr = std::make_shared<ode_data_config>(ode_coeffs_data_ptr, ode_nonhom_data_ptr);
    // boundary conditions:
    auto const &dirichlet = 1.0;
    auto const &robin_first = 1.0;
    auto const &robin_second = 0.0;
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_high_ptr = std::make_shared<robin_boundary_1d>(robin_first, robin_second);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize ode solver
    ode_equation odesolver(ode_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                           dev_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    odesolver.solve(solution);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + 1.0); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-4, "element at " << t << "must be appr same.");
    }
}

void test_impl_simple_ode_dirichlet_robin_bc_cuda_solver_device()
{
    std::cout << "==================================================\n";
    std::cout << "================= (Dir-Rob BC) ===================\n";
    std::cout << "==================================================\n";

    impl_simple_ode_dirichlet_robin_bc_cuda_solver_device_qr();

    std::cout << "==================================================\n";
}

void impl_simple_ode_neumann_robin_bc_cuda_solver_device_qr()
{

    std::cout << "=================================\n";
    std::cout << "Solving Boundary-value problem: \n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " u''(t) = 6*t, \n\n";
    std::cout << " where\n\n";
    std::cout << " t in <0,2>,\n";
    std::cout << " u'(0) = 0 \n";
    std::cout << " u'(2) + 2*u(t) = 0\n\n";
    std::cout << "Exact solution is:\n\n";
    std::cout << " u(t) = t*t*t - 14\n";
    std::cout << "=================================\n";

    // number of space subdivisions:
    std::size_t Sd{100};
    // space range:
    auto const space_range = std::make_shared<range>(0.0, 2.0);
    // discretization config:
    auto const discretization_ptr = std::make_shared<ode_discretization_config>(space_range, Sd);
    // coeffs:
    auto a = [](double x) { return 0.0; };
    auto b = [](double x) { return 0.0; };
    auto const ode_coeffs_data_ptr = std::make_shared<ode_coefficient_data_config>(a, b);
    // nonhom data:
    auto two = [](double x) { return 6 * x; };
    auto const ode_nonhom_data_ptr = std::make_shared<ode_nonhom_data_config>(two);
    // ode data config:
    auto const ode_data_ptr = std::make_shared<ode_data_config>(ode_coeffs_data_ptr, ode_nonhom_data_ptr);
    // boundary conditions:
    auto const &neumann = 0.0;
    auto const &robin_first = 2.0;
    auto const &robin_second = 0.0;
    auto const &boundary_low_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_high_ptr = std::make_shared<robin_boundary_1d>(robin_first, robin_second);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize ode solver
    ode_equation odesolver(ode_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                           dev_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    odesolver.solve(solution);

    // exact value:
    auto exact = [](double x) { return (x * x * x - static_cast<double>(14.0)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-2, "element at " << t << "must be appr same.");
    }
}

void test_impl_simple_ode_neumann_robin_bc_cuda_solver_device()
{
    std::cout << "==================================================\n";
    std::cout << "=================  (Neu-Rob BC) ==================\n";
    std::cout << "==================================================\n";

    impl_simple_ode_neumann_robin_bc_cuda_solver_device_qr();

    std::cout << "==================================================\n";
}

void impl_simple_ode1_neumann_robin_bc_cuda_solver_device_qr()
{

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

    // number of space subdivisions:
    std::size_t Sd{100};
    // space range:
    auto const space_range = std::make_shared<range>(0.0, 1.0);
    // discretization config:
    auto const discretization_ptr = std::make_shared<ode_discretization_config>(space_range, Sd);
    // coeffs:
    auto a = [](double x) { return 0.0; };
    auto b = [](double x) { return 0.0; };
    auto const ode_coeffs_data_ptr = std::make_shared<ode_coefficient_data_config>(a, b);
    // nonhom data:
    auto two = [](double x) { return -2.0; };
    auto const ode_nonhom_data_ptr = std::make_shared<ode_nonhom_data_config>(two);
    // ode data config:
    auto const ode_data_ptr = std::make_shared<ode_data_config>(ode_coeffs_data_ptr, ode_nonhom_data_ptr);
    // boundary conditions:
    auto const &neumann = -1.0;
    auto const &robin_first = 2.0;
    auto const &robin_second = 0.0;
    auto const &boundary_low_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_high_ptr = std::make_shared<robin_boundary_1d>(robin_first, robin_second);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize ode solver
    ode_equation odesolver(ode_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                           dev_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    odesolver.solve(solution);

    // exact value:
    auto exact = [](double x) { return (-x * x + x + (0.5)); };

    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element at " << t << "must be appr same.");
    }
}

void test_impl_simple_ode1_neumann_robin_bc_cuda_solver_device()
{
    std::cout << "==================================================\n";
    std::cout << "================ (Neu - Rob) =====================\n";
    std::cout << "==================================================\n";

    impl_simple_ode1_neumann_robin_bc_cuda_solver_device_qr();

    std::cout << "==================================================\n";
}
