#include "lss_advection_equation_t.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../discretization/lss_grid.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../../discretization/lss_grid_config_hints.hpp"
#include "../../../pde_solvers/1d/heat_type/lss_heat_equation.hpp"
#include "../../../pde_solvers/lss_heat_data_config.hpp"
#include "../../../pde_solvers/lss_heat_solver_config.hpp"
#include <map>

using lss_boundary::dirichlet_boundary_1d;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_grids::grid_1d;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d;
using lss_grids::grid_transform_config_1d;
using lss_pde_solvers::heat_coefficient_data_config_1d;
using lss_pde_solvers::heat_data_config_1d;
using lss_pde_solvers::heat_implicit_solver_config;
using lss_pde_solvers::heat_initial_data_config_1d;
using lss_pde_solvers::pde_discretization_config_1d;
using lss_pde_solvers::one_dimensional::implicit_solvers::heat_equation;
using lss_utility::container_t;
using lss_utility::pi;
using lss_utility::range;

// ///////////////////////////////////////////////////////////////////////////
//							ADVECTION PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// ====== Advection Diffusion problem with homogeneous boundary conditions ===
// ===========================================================================

void impl_ddv_diff_equation_dirichlet_bc_cuda_solver_device_qr_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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

void impl_adv_diff_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) and implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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

void test_impl_adv_diff_equation_dirichlet_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Advection (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_ddv_diff_equation_dirichlet_bc_cuda_solver_device_qr_euler();
    impl_adv_diff_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_adv_diff_equation_dirichlet_bc_sor_solver_device_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using CUDA SOR solver (DEVICE) with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
    details["sor_omega"] = 1.0;
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_adv_diff_equation_dirichlet_bc_sor_solver_device_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using CUDA SOR solver (DEVICE) with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
    details["sor_omega"] = 1.0;
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_adv_diff_equation_dirichlet_bc_sor_solver_device()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Advection  (SOR QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_adv_diff_equation_dirichlet_bc_sor_solver_device_euler();
    impl_adv_diff_equation_dirichlet_bc_sor_solver_device_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_adv_diff_equation_dirichlet_bc_sor_solver_host_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using SOR solver with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
    details["sor_omega"] = 1.0;
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_adv_diff_equation_dirichlet_bc_sor_solver_host_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using SOR solver with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
    details["sor_omega"] = 1.0;
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_adv_diff_equation_dirichlet_bc_sor_solver_host()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Advection (SOR QR HOST) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    impl_adv_diff_equation_dirichlet_bc_sor_solver_host_euler();
    impl_adv_diff_equation_dirichlet_bc_sor_solver_host_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using CUDA Solver on HOST with QR and implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using CUDA Solver on HOST with QR and implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Advection(CUDA QR HOST) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr_euler();
    impl_adv_diff_equation_dirichlet_bc_cuda_solver_host_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_adv_diff_equation_dirichlet_bc_double_sweep_solver_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_dssolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using Double Sweep on HOST with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_adv_diff_equation_dirichlet_bc_double_sweep_solver_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_fwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using Double Sweep on HOST with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_adv_diff_equation_dirichlet_bc_double_sweep_solver()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Advection(Double Sweep) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    impl_adv_diff_equation_dirichlet_bc_double_sweep_solver_euler();
    impl_adv_diff_equation_dirichlet_bc_double_sweep_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_tlusolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using Thomas LU algorithm with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::host_fwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Advection Diffusion equation: \n\n";
    std::cout << " Using Thomas LU algorithm with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t)  +  U_x(x,t) = U_xx(x,t), \n\n";
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
    auto b = [](double t, double x) { return -1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, other);
    // initial condition:
    auto initial_condition = [](double x) { return 1.0; };
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
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const first = 2.0 / pi();
        double const exp_0p5x = std::exp(0.5 * x);
        double const exp_m0p5 = std::exp(-0.5);
        double np_sqr{};
        double sum{};
        double num{}, den{}, var{};
        double lambda{};
        for (std::size_t i = 1; i <= n; ++i)
        {
            np_sqr = (i * i * pi() * pi());
            lambda = 0.25 + np_sqr;
            num =
                (1.0 - std::pow(-1.0, i) * exp_m0p5) * exp_0p5x * std::exp(-1.0 * lambda * t) * std::sin(i * pi() * x);
            den = i * (1.0 + (0.25 / np_sqr));
            var = num / den;
            sum += var;
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
        benchmark = exact(x, time_range->upper(), 20);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Advection (Thomas LU) Equation (Dirichlet BC) ==\n";
    std::cout << "============================================================\n";

    impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver_euler();
    impl_adv_diff_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson();

    std::cout << "============================================================\n";
}
