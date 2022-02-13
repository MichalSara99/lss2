#include "lss_pure_heat_equation_t.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../discretization/lss_grid.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../../discretization/lss_grid_config_hints.hpp"
#include "../../../pde_solvers/1d/heat_type/lss_heat_equation.hpp"
#include "../../../pde_solvers/lss_heat_data_config.hpp"
#include "../../../pde_solvers/lss_heat_solver_config.hpp"
#include <map>

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using lss_containers::matrix_2d;
using lss_enumerations::grid_enum;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_grids::grid_1d;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d;
using lss_grids::grid_transform_config_1d;
using lss_pde_solvers::heat_coefficient_data_config_1d;
using lss_pde_solvers::heat_data_config_1d;
using lss_pde_solvers::heat_implicit_solver_config;
using lss_pde_solvers::heat_initial_data_config_1d;
using lss_pde_solvers::heat_source_data_config_1d;
using lss_pde_solvers::pde_discretization_config_1d;
using lss_pde_solvers::one_dimensional::implicit_solvers::heat_equation;
using lss_utility::container_t;
using lss_utility::pi;
using lss_utility::range;

// ///////////////////////////////////////////////////////////////////////////
//							PURE HEAT PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Heat problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = 2.0 / pi();
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }

    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Heat (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_euler();
    impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_dirichlet_bc_sor_solver_device_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA SOR solver (DEVICE) with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_dirichlet_bc_sor_solver_device_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA SOR solver (DEVICE) with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_sor_solver_device()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Heat (SOR QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_sor_solver_device_euler();
    impl_pure_heat_equation_dirichlet_bc_sor_solver_device_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_dirichlet_bc_sor_solver_host_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using SOR solver with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_dirichlet_bc_sor_solver_host_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using SOR solver with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_sor_solver_host()
{

    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Heat (SOR QR HOST) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_sor_solver_host_euler();
    impl_pure_heat_equation_dirichlet_bc_sor_solver_host_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA Solver on HOST with QR and implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA Solver on HOST with QR and implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time rang
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr()
{

    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Heat(CUDA QR HOST) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr_euler();
    impl_pure_heat_equation_dirichlet_bc_cuda_solver_host_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_dirichlet_bc_double_sweep_solver_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_dssolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Double Sweep on HOST with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_dirichlet_bc_double_sweep_solver_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Double Sweep on HOST with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_double_sweep_solver()
{

    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Heat(Double Sweep) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_double_sweep_solver_euler();
    impl_pure_heat_equation_dirichlet_bc_double_sweep_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_tlusolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Thomas LU algorithm with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Thomas LU algorithm with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    // note: size is Sd+1 since we must include space point at x = 0
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver()
{

    std::cout << "============================================================\n";
    std::cout << "== Implicit Pure Heat (Thomas LU) Equation (Dirichlet BC) ==\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver_euler();
    impl_pure_heat_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

// Neuman-Dirichlet Boundaries:

void impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA Solver on DEVICE with QR with implicit Euler\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = 0 ,U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = cos((pi()/2)*x), x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::cos(x * pi() * static_cast<double>(0.5)); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_left_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_right_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_left_ptr, boundary_right_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    // get exact solution:
    auto exact = [](double x, double t) {
        double const pipi = (pi() * pi());
        double const expon = (-pipi * 0.25);
        double res = std::exp(expon * t) * std::cos(x * pi() * 0.5);
        return res;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper());
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper());
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA Solver on DEVICE with QR with implicit CN\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = 0 ,U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = cos((pi()/2)*x), x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::cos(x * pi() * static_cast<double>(0.5)); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_lower_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_upper_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_lower_ptr, boundary_upper_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:

    auto exact = [](double x, double t) {
        double const pipi = (pi() * pi());
        double const expon = (-pipi * 0.25);
        double res = std::exp(expon * t) * std::cos(x * pi() * 0.5);
        return res;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper());
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper());
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Pure Heat (CUDA DEV QR) Equation (Neu-Dir BC) ==\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr_euler();
    impl_pure_heat_equation_neumann_bc_cuda_solver_device_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

// Neumann-Neumann Boundaries:
void impl_pure_heat_equation_neumann_bc_thomas_lu_solver_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_tlusolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Thomas LU solver with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = U_x(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const pipi = (pi() * pi());
        double const first = (4.0) / pipi;
        double sum{};
        double var0{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var0 = static_cast<double>(2 * i - 1);
            var1 = std::exp(-1.0 * pipi * var0 * var0 * t);
            var2 = std::cos(var0 * pi() * x) / (var0 * var0);
            sum += (var1 * var2);
        }
        return (0.5 - first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}
void impl_pure_heat_equation_neumann_bc_thomas_lu_solver_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Thomas LU solver with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = U_x(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.3);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const pipi = (pi() * pi());
        double const first = (4.0) / pipi;
        double sum{};
        double var0{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var0 = static_cast<double>(2 * i - 1);
            var1 = std::exp(-1.0 * pipi * var0 * var0 * t);
            var2 = std::cos(var0 * pi() * x) / (var0 * var0);
            sum += (var1 * var2);
        }
        return (0.5 - first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_neumann_bc_thomas_lu_solver()
{
    std::cout << "============================================================\n";
    std::cout << "=== Implicit Pure Heat (Thomas LU) Equation (Neu-Neu BC) ===\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_neumann_bc_thomas_lu_solver_euler();
    impl_pure_heat_equation_neumann_bc_thomas_lu_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_neumann_bc_double_sweep_solver_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_dssolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Double Sweep algorithm with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = U_x(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const pipi = (pi() * pi());
        double const first = (4.0) / pipi;
        double sum{};
        double var0{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var0 = static_cast<double>(2 * i - 1);
            var1 = std::exp((-1.0) * pipi * var0 * var0 * t);
            var2 = std::cos(var0 * pi() * x) / (var0 * var0);
            sum += (var1 * var2);
        }
        return ((0.5) - first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_neumann_bc_double_sweep_solver_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Double Sweep algorithm with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = U_x(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const pipi = (pi() * pi());
        double const first = (4.0) / pipi;
        double sum{};
        double var0{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var0 = static_cast<double>(2 * i - 1);
            var1 = std::exp((-1.0) * pipi * var0 * var0 * t);
            var2 = std::cos(var0 * pi() * x) / (var0 * var0);
            sum += (var1 * var2);
        }
        return ((0.5) - first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_neumann_bc_double_sweep_solver()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Heat (Double Sweep) Equation (Neu-Neu BC) ==\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_neumann_bc_double_sweep_solver_euler();
    impl_pure_heat_equation_neumann_bc_double_sweep_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

// get the whole surface with stepping:
void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_euler_stepping()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double const k = discretization_ptr->time_step();
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        std::cout << "time: " << t * k << ":\n";
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark = exact(x, t * k, 20);
            std::cout << "t_" << j << ": " << solutions(t, j) << " |  " << benchmark << " | "
                      << (solutions(t, j) - benchmark) << '\n';
        }
    }
}

void impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson_stepping()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double const k = discretization_ptr->time_step();
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        std::cout << "time: " << t * k << ":\n";
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark = exact(x, t * k, 20);
            std::cout << "t_" << j << ": " << solutions(t, j) << " |  " << benchmark << " | "
                      << (solutions(t, j) - benchmark) << '\n';
        }
    }
}

void test_impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_stepping()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Heat (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_euler_stepping();
    impl_pure_heat_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson_stepping();

    std::cout << "============================================================\n";
}

// ===========================================================================
// ====== Heat problem with homogeneous boundary conditions & source =========
// ===========================================================================

void impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device_qr_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t) + x, \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = 1, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat source data:
    auto source = [](double t, double x) { return x; };
    auto const heat_source_data_ptr = std::make_shared<heat_source_data_config_1d>(source);
    // heat data config:
    auto const heat_data_ptr =
        std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr, heat_source_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);

    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double sum{};
        double q_n{};
        double f_n{};
        double lam_n{};
        double lam_2{};
        double var1{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            q_n = (2.0 / (i * pi())) * std::pow(-1.0, i + 1);
            f_n = (2.0 / (i * pi())) * (1.0 - std::pow(-1.0, i));
            lam_n = i * pi();
            lam_2 = lam_n * lam_n;
            var1 = ((q_n / lam_2) + (f_n - (q_n / lam_2)) * std::exp(-1.0 * lam_2 * t)) * std::sin(i * pi() * x);
            sum += var1;
        }
        return sum;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 30);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device_qr_clark_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t) + x, \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = 1, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat source data:
    auto source = [](double t, double x) { return x; };
    auto const heat_source_data_ptr = std::make_shared<heat_source_data_config_1d>(source);
    // heat data config:
    auto const heat_data_ptr =
        std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr, heat_source_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);

    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double sum{};
        double q_n{};
        double f_n{};
        double lam_n{};
        double lam_2{};
        double var1{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            q_n = (2.0 / (i * pi())) * std::pow(-1.0, i + 1);
            f_n = (2.0 / (i * pi())) * (1.0 - std::pow(-1.0, i));
            lam_n = i * pi();
            lam_2 = lam_n * lam_n;
            var1 = ((q_n / lam_2) + (f_n - (q_n / lam_2)) * std::exp(-1.0 * lam_2 * t)) * std::sin(i * pi() * x);
            sum += var1;
        }
        return sum;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 30);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Heat (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device_qr_euler();
    impl_pure_heat_equation_source_dirichlet_bc_cuda_solver_device_qr_clark_nicolson();

    std::cout << "============================================================\n";
}

void impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA SOR solver (DEVICE) with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t) + x, \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = 1, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat source data:
    auto source = [](double t, double x) { return x; };
    auto const heat_source_data_ptr = std::make_shared<heat_source_data_config_1d>(source);
    // heat data config:
    auto const heat_data_ptr =
        std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr, heat_source_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);

    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double sum{};
        double q_n{};
        double f_n{};
        double lam_n{};
        double lam_2{};
        double var1{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            q_n = (2.0 / (i * pi())) * std::pow(-1.0, i + 1);
            f_n = (2.0 / (i * pi())) * (1.0 - std::pow(-1.0, i));
            lam_n = i * pi();
            lam_2 = lam_n * lam_n;
            var1 = ((q_n / lam_2) + (f_n - (q_n / lam_2)) * std::exp(-1.0 * lam_2 * t)) * std::sin(i * pi() * x);
            sum += var1;
        }
        return sum;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 30);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device_clark_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using CUDA SOR solver (DEVICE) with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t) + x, \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = 1, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat source data:
    auto source = [](double t, double x) { return x; };
    auto const heat_source_data_ptr = std::make_shared<heat_source_data_config_1d>(source);
    // heat data config:
    auto const heat_data_ptr =
        std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr, heat_source_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);

    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double sum{};
        double q_n{};
        double f_n{};
        double lam_n{};
        double lam_2{};
        double var1{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            q_n = (2.0 / (i * pi())) * std::pow(-1.0, i + 1);
            f_n = (2.0 / (i * pi())) * (1.0 - std::pow(-1.0, i));
            lam_n = i * pi();
            lam_2 = lam_n * lam_n;
            var1 = ((q_n / lam_2) + (f_n - (q_n / lam_2)) * std::exp(-1.0 * lam_2 * t)) * std::sin(i * pi() * x);
            sum += var1;
        }
        return sum;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 30);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Heat (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device_euler();
    impl_pure_heat_equation_source_dirichlet_bc_sor_solver_device_clark_nicolson();

    std::cout << "============================================================\n";
}

// ===========================================================================
// ========================== EXPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Heat problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

void expl_pure_heat_equation_dirichlet_bc_barakat_clark()
{

    using lss_pde_solvers::default_heat_solver_configs::host_expl_fwd_bc_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using explicit Barakat-Clark method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 300;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_bc_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void expl_pure_heat_equation_dirichlet_bc_saulyev()
{
    using lss_pde_solvers::default_heat_solver_configs::host_expl_fwd_s_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;
    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using explicit Saulyev method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 300;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_s_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void expl_pure_heat_equation_dirichlet_bc_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_expl_fwd_euler_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using explicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 10000;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = static_cast<double>(2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_expl_pure_heat_equation_dirichlet_bc_ade()
{
    std::cout << "============================================================\n";
    std::cout << "===== Explicit Pure Heat (ADE ) Equation (Dirichlet BC) ====\n";
    std::cout << "============================================================\n";

    expl_pure_heat_equation_dirichlet_bc_barakat_clark();
    expl_pure_heat_equation_dirichlet_bc_saulyev();
    expl_pure_heat_equation_dirichlet_bc_euler();

    std::cout << "============================================================\n";
}

void empl_pure_heat_equation_neumann_dirichlet_bc_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_expl_fwd_euler_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = 0 ,U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = cos((pi()/2)*x), x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 10000;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::cos(x * pi() * static_cast<double>(0.5)); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_lower_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_upper_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_lower_ptr, boundary_upper_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t) {
        double const pipi = (pi() * pi());
        double const expon = (-pipi * (0.25));
        double res = std::exp(expon * t) * std::cos(x * pi() * (0.5));
        return res;
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper());
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper());
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void expl_pure_heat_equation_neumann_neumann_bc_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_expl_fwd_euler_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U_x(0,t) = U_x(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 10000;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.3);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const pipi = (pi() * pi());
        double const first = (4.0) / pipi;
        double sum{};
        double var0{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var0 = static_cast<double>(2 * i - 1);
            var1 = std::exp((-1.0) * pipi * var0 * var0 * t);
            var2 = std::cos(var0 * pi() * x) / (var0 * var0);
            sum += (var1 * var2);
        }
        return ((0.5) - first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_expl_pure_heat_equation_neumann_bc_euler()
{
    std::cout << "============================================================\n";
    std::cout << "===== Explicit Pure Heat (ADE ) Equation (Non-Dir BC) ======\n";
    std::cout << "============================================================\n";

    empl_pure_heat_equation_neumann_dirichlet_bc_euler();
    expl_pure_heat_equation_neumann_neumann_bc_euler();

    std::cout << "============================================================\n";
}

void expl_pure_heat_equation_dirichlet_bc_euler_device()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_expl_fwd_euler_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heat equation: \n\n";
    std::cout << " Using explicit Euler (on DEV) method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 10000;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.1);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return x; };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(initial_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_expl_fwd_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = (2.0 / pi());
        double sum{};
        double var1{};
        double var2{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            var1 = std::pow(-1.0, i + 1) * std::exp(-1.0 * (i * pi()) * (i * pi()) * t);
            var2 = std::sin(i * pi() * x) / i;
            sum += (var1 * var2);
        }
        return (first * sum);
    };

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = exact(x, time_range->upper(), 20);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = exact(x, time_range->upper(), 30);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-3, "element " << t << " must be approx same");
    }
}

void test_expl_pure_heat_equation_dirichlet_bc_device()
{
    std::cout << "============================================================\n";
    std::cout << "===== Explicit Pure Heat (Device ) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    expl_pure_heat_equation_dirichlet_bc_euler_device();

    std::cout << "============================================================\n";
}
