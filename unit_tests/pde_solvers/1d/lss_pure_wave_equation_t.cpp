#include "lss_pure_wave_equation_t.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../pde_solvers/1d/wave_type/lss_wave_equation.hpp"
#include <map>

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using lss_enumerations::grid_enum;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_grids::grid_1d;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d;
using lss_grids::grid_transform_config_1d;
using lss_pde_solvers::pde_discretization_config_1d;
using lss_pde_solvers::wave_coefficient_data_config_1d;
using lss_pde_solvers::wave_data_config_1d;
using lss_pde_solvers::wave_implicit_solver_config;
using lss_pde_solvers::wave_initial_data_config_1d;
using lss_pde_solvers::default_wave_solver_configs::dev_fwd_cusolver_qr_solver_config_ptr;
using lss_pde_solvers::one_dimensional::implicit_solvers::wave_equation;
using lss_utility::container_t;
using lss_utility::pi;
using lss_utility::range;

// ///////////////////////////////////////////////////////////////////////////
//							PURE WAVE PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Wave problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

void impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_qr_detail()
{
    using lss_pde_solvers::default_wave_solver_configs::dev_fwd_cusolver_qr_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 0.8);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>((0.5), alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-4, "element " << t << " must be approx same");
    }
}

void test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Wave (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_qr_detail();

    std::cout << "============================================================\n";
}

void impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_sor_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::dev_fwd_sorsolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using CUDA solver with LU (DEVICE) method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.7));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = 1.0;
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>((0.5), alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_sorsolver_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-4, "element " << t << " must be approx same");
    }
}

void test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_sor()
{
    std::cout << "============================================================\n";
    std::cout << "Implicit Pure Wave (CUDA SOR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_wave_equation_dirichlet_bc_cuda_solver_device_sor_detail();

    std::cout << "============================================================\n";
}

void impl_pure_wave_equation_dirichlet_bc_cuda_solver_host_sor_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::host_fwd_sorsolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using CUDA solver with LU (DEVICE) method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.7));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = 1.0;
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_sorsolver_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-4, "element " << t << " must be approx same");
    }
}

void test_impl_pure_wave_equation_dirichlet_bc_cuda_solver_host_sor()
{
    std::cout << "============================================================\n";
    std::cout << "Implicit Pure Wave (CUDA SOR Host) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_wave_equation_dirichlet_bc_cuda_solver_host_sor_detail();

    std::cout << "============================================================\n";
}

void impl_pure_wave_equation_dirichlet_bc_solver_host_double_sweep_detail()
{
    using lss_pde_solvers::default_wave_solver_configs::host_fwd_dssolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using Double Sweep Solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.7));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-4, "element " << t << " must be approx same");
    }
}

void test_impl_pure_wave_equation_dirichlet_bc_solver_host_double_sweep()
{
    std::cout << "============================================================\n";
    std::cout << "Implicit Pure Wave (Double Sweep) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_pure_wave_equation_dirichlet_bc_solver_host_double_sweep_detail();

    std::cout << "============================================================\n";
}

void impl_pure_wave_equation_dirichlet_bc_solver_host_lu_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::host_fwd_tlusolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using Thomas LU Solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.7));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-4, "element " << t << " must be approx same");
    }
}

void test_impl_pure_wave_equation_dirichlet_bc_solver_host_lu()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Pure Wave (Thomas LU) Equation (Dirichlet BC) ==\n";
    std::cout << "============================================================\n";

    impl_pure_wave_equation_dirichlet_bc_solver_host_lu_detail();

    std::cout << "============================================================\n";
}

void impl_wave_equation_dirichlet_bc_solver_host_lu_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::host_fwd_tlusolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using Thomas LU Solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = 4*U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = -4*t*t, t > 0 \n\n";
    std::cout << " U(x,0) = x*(1 - x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 8*x, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.7));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 4.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto first_initial_condition = [](double x) { return (x * (1.0 - x)); };
    auto second_initial_condition = [](double x) { return (8 * x); };
    auto const wave_init_data_ptr =
        std::make_shared<wave_initial_data_config_1d>(first_initial_condition, second_initial_condition);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_0 = [](double t) { return (-4.0 * t * t); };
    auto const &dirichlet_1 = [](double t) { return (-4.0 * t * t + 8.0 * t); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_0);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_1);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double res = x - x * x - 4.0 * t * t + 8.0 * t * x;
        return (res);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 2e-4, "element " << t << " must be approx same");
    }
}

void test_impl_wave_equation_dirichlet_bc_solver_host_lu()
{
    std::cout << "============================================================\n";
    std::cout << "===== Implicit Wave (Thomas LU) Equation (Dirichlet BC) ====\n";
    std::cout << "============================================================\n";

    impl_wave_equation_dirichlet_bc_solver_host_lu_detail();

    std::cout << "============================================================\n";
}

void impl_damped_wave_equation_dirichlet_bc_solver_host_double_sweep_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::host_fwd_dssolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using Double Sweep Solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t)  + U_t(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,pi> and t > 0,\n";
    std::cout << " U(0,t) = U(pi,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(pi()));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.7));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [](double t, double x) { return 1.0; };
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(a, b, other, other);
    // initial condition:
    auto first_initial_condition = [](double x) { return std::sin(x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(first_initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_0 = [](double t) { return 0.0; };
    auto const &dirichlet_1 = [](double t) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_0);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_1);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double exp_half = std::exp(-0.5 * t);
        const double sqrt_3 = std::sqrt(3.0);
        const double arg = 0.5 * sqrt_3;
        const double res = exp_half * std::sin(x) * (std::cos(arg * t) + (sin(arg * t) / sqrt_3));
        return (res);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-4, "element " << t << " must be approx same");
    }
}

void test_impl_damped_wave_equation_dirichlet_bc_solver_host_double_sweep()
{
    std::cout << "============================================================\n";
    std::cout << "==== Implicit Wave (Double Sweep) Equation (Dirichlet BC) ==\n";
    std::cout << "============================================================\n";

    impl_damped_wave_equation_dirichlet_bc_solver_host_double_sweep_detail();

    std::cout << "============================================================\n";
}

// Neumann BC:

void impl_pure_wave_equation_neumann_bc_cuda_solver_device_qr_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::dev_fwd_cusolver_qr_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using CUDA solver with QR (DEVICE) method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = 4*U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,pi> and t > 0,\n";
    std::cout << " U_x(0,t) = U_x(pi,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = 3*cos(x), x in <0,pi> \n\n";
    std::cout << " U_x(x,0) = 1 - cos(4*x), x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(pi()));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.8));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 4.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto first_initial_condition = [](double x) { return 3.0 * std::cos(x); };
    auto second_initial_condition = [](double x) { return (1.0 - std::cos(4.0 * x)); };
    auto const wave_init_data_ptr =
        std::make_shared<wave_initial_data_config_1d>(first_initial_condition, second_initial_condition);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &neumann = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<neumann_boundary_1d>(neumann);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double var1 = 3.0 * std::cos(2.0 * t) * std::cos(x);
        const double var2 = -0.125 * std::sin(8.0 * t) * std::cos(4.0 * x);
        return (t + var1 + var2);
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
}

void test_impl_pure_wave_equation_neumann_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Wave (CUDA QR DEVICE) Equation (Neumann BC)=\n";
    std::cout << "============================================================\n";

    impl_pure_wave_equation_neumann_bc_cuda_solver_device_qr_detail();

    std::cout << "============================================================\n";
}

// ===========================================================================
// ========================== EXPLICIT SOLVERS ===============================
// ===========================================================================

// ===========================================================================
// =========== Wave problem with homogeneous boundary conditions =============
// ===========================================================================

// Dirichlet boundaries:

void expl_pure_wave_equation_dirichlet_bc_cuda_host_solver_detail()
{
    using lss_pde_solvers::default_wave_solver_configs::host_expl_fwd_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::wave_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using CUDA solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 50;
    // number of time subdivisions:
    std::size_t const Td = 2000;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.8));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-4, "element " << t << " must be approx same");
    }
}

void test_expl_pure_wave_equation_dirichlet_bc_cuda_host_solver()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Wave (CUDA DEVICE) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    expl_pure_wave_equation_dirichlet_bc_cuda_host_solver_detail();

    std::cout << "============================================================\n";
}

void expl_pure_wave_equation_dirichlet_bc_cuda_device_solver_detail()
{

    using lss_pde_solvers::default_wave_solver_configs::dev_expl_fwd_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::wave_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using CUDA solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = U(1,t) = 0, t > 0 \n\n";
    std::cout << " U(x,0) = sin(pi*x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 50;
    // number of time subdivisions:
    std::size_t const Td = 2000;
    // space range:
    auto const &space_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(1.0));
    // time range
    auto const &time_range = std::make_shared<range>(static_cast<double>(0.0), static_cast<double>(0.8));
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto b = [](double t, double x) { return 1.0; };
    auto other = [](double t, double x) { return 0.0; };
    auto const wave_coeffs_data_ptr = std::make_shared<wave_coefficient_data_config_1d>(other, b, other, other);
    // initial condition:
    auto initial_condition = [](double x) { return std::sin(pi() * x); };
    auto zero = [](double x) { return 0.0; };
    auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, zero);
    // wave data config:
    auto const wave_data_ptr = std::make_shared<wave_data_config_1d>(wave_coeffs_data_ptr, wave_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet = [](double t) { return 0.0; };
    auto const &boundary_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet);
    auto const &boundary_pair = std::make_pair(boundary_ptr, boundary_ptr);
    // grid:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_expl_fwd_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
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
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-4, "element " << t << " must be approx same");
    }
}

void test_expl_pure_wave_equation_dirichlet_bc_cuda_device_solver()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Wave (CUDA DEVICE) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    expl_pure_wave_equation_dirichlet_bc_cuda_device_solver_detail();

    std::cout << "============================================================\n";
}
