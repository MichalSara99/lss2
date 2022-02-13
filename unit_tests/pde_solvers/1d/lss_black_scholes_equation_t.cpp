#include "lss_black_scholes_equation_t.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
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
using lss_pde_solvers::pde_discretization_config_1d;
using lss_pde_solvers::one_dimensional::implicit_solvers::heat_equation;
using lss_utility::black_scholes_exact;
using lss_utility::container_t;
using lss_utility::range;

// ///////////////////////////////////////////////////////////////////////////
//							BLACK SCHOLES PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

// Dirichlet boundaries:
void impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using CUDA on DEVICE with QR implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_bwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using CUDA on DEVICE with QR implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_bwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Black-Scholes (CUDA QR DEVICE) Equation (Dir BC) \n";
    std::cout << "============================================================\n";

    impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_euler();
    impl_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_black_scholes_equation_dirichlet_bc_sor_solver_device_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using SOR on DEVICE with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_bwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_black_scholes_equation_dirichlet_bc_sor_solver_device_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using SOR on DEVICE with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_bwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_black_scholes_equation_dirichlet_bc_sor_solver_device()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Black-Scholes (SOR DEVICE) Equation (Dir BC) ====\n";
    std::cout << "============================================================\n";

    impl_black_scholes_equation_dirichlet_bc_sor_solver_device_euler();
    impl_black_scholes_equation_dirichlet_bc_sor_solver_device_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_black_scholes_equation_dirichlet_bc_sor_solver_host_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_sorsolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using SOR on HOST with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_sorsolver_euler_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_black_scholes_equation_dirichlet_bc_sor_solver_host_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_sorsolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using SOR on HOST with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // details:
    std::map<std::string, double> details;
    details["sor_omega"] = static_cast<double>(1.0);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_sorsolver_cn_solver_config_ptr, details);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_black_scholes_equation_dirichlet_bc_sor_solver_host()
{
    std::cout << "============================================================\n";
    std::cout << "=== Implicit Black-Scholes (SOR HOST) Equation (Dir BC) ====\n";
    std::cout << "============================================================\n";

    impl_black_scholes_equation_dirichlet_bc_sor_solver_host_euler();
    impl_black_scholes_equation_dirichlet_bc_sor_solver_host_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_black_scholes_equation_dirichlet_bc_double_sweep_solver_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Double Sweep on HOST with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_dssolver_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_black_scholes_equation_dirichlet_bc_double_sweep_solver_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Double Sweep on HOST with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const alpha_scale = 3.0;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_black_scholes_equation_dirichlet_bc_double_sweep_solver()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Black-Scholes (Double Sweep) Equation (Dir BC) ==\n";
    std::cout << "============================================================\n";

    impl_black_scholes_equation_dirichlet_bc_double_sweep_solver_euler();
    impl_black_scholes_equation_dirichlet_bc_double_sweep_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Thomas LU on HOST with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_tlusolver_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-3, "element " << t << " must be approx same");
    }
}

void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Thomas LU on HOST with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-3, "element " << t << " must be approx same");
    }
}

void test_impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Black-Scholes (Thomas LU) Equation (Dir BC) ====\n";
    std::cout << "============================================================\n";

    impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_euler();
    impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

// forward starting call:
void impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_euler()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_cusolver_qr_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using CUDA on DEVICE with QR implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0.5 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0.5 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    auto const &fwd_start = 0.5;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(fwd_start, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_bwd_cusolver_qr_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(fwd_start, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using CUDA on DEVICE with QR implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0.5 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0.5 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    auto const &fwd_start = 0.5;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(fwd_start, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_bwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(fwd_start, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Black-Scholes (CUDA QR DEVICE) Equation (Dir BC) \n";
    std::cout << "============================================================\n";

    impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_euler();
    impl_fwd_black_scholes_equation_dirichlet_bc_cuda_solver_device_qr_crank_nicolson();

    std::cout << "============================================================\n";
}

// Uisng stepping = getting the whole surface
void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_euler_stepping()
{

    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_euler_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Thomas LU on HOST with implicit Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_tlusolver_euler_solver_config_ptr);
    // prepare container for solutions:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double const k = discretization_ptr->time_step();
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    // for (std::size_t t = 0; t < solutions.rows(); ++t)
    //{
    //     std::cout << "time: " << t * k << ":\n";
    //     for (std::size_t j = 0; j < solutions.columns(); ++j)
    //     {
    //         x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
    //         benchmark = bs_exact.call(x, maturity - t * k);
    //         std::cout << "t_" << j << ": " << solutions(t, j) << " |  " << benchmark << " | "
    //                   << (solutions(t, j) - benchmark) << '\n';
    //     }
    // }
    //  TEST:
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        std::cout << "time: " << t * k << ":\n";
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark = bs_exact.call(x, maturity - t * k);
            LSS_ASSERT(std::abs(solutions(t, j) - benchmark) < 3e-2,
                       "element (" << t << "," << j << ") must be approx same");
        }
    }
}

void impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson_stepping()
{

    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Thomas LU on HOST with implicit CN method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solutions:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double const k = discretization_ptr->time_step();
    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    // for (std::size_t t = 0; t < solutions.rows(); ++t)
    //{
    //     std::cout << "time: " << t * k << ":\n";
    //     for (std::size_t j = 0; j < solutions.columns(); ++j)
    //     {
    //         x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
    //         benchmark = bs_exact.call(x, maturity - t * k);
    //         std::cout << "t_" << j << ": " << solutions(t, j) << " |  " << benchmark << " | "
    //                   << (solutions(t, j) - benchmark) << '\n';
    //     }
    // }

    //  TEST:
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        std::cout << "time: " << t * k << ":\n";
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark = bs_exact.call(x, maturity - t * k);
            LSS_ASSERT(std::abs(solutions(t, j) - benchmark) < 8e-3,
                       "element (" << t << "," << j << ") must be approx same");
        }
    }
}

void test_impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_stepping()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Black-Scholes (Thomas LU) Equation (Dir BC) ====\n";
    std::cout << "============================================================\n";

    impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_euler_stepping();
    impl_black_scholes_equation_dirichlet_bc_thomas_lu_solver_crank_nicolson_stepping();

    std::cout << "============================================================\n";
}

// ===========================================================================
// ========================== EXPLICIT SOLVERS ===============================
// ===========================================================================

// Dirichlet boundaries:

void expl_black_scholes_equation_dirichlet_bc_barakat_clark()
{
    using lss_pde_solvers::default_heat_solver_configs::host_expl_bwd_bc_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Barakat-Clark ADE method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 300;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_bwd_bc_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 4e-3, "element " << t << " must be approx same");
    }
}

void expl_black_scholes_equation_dirichlet_bc_saulyev()
{

    using lss_pde_solvers::default_heat_solver_configs::host_expl_bwd_s_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Saulyev ADE method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 300;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_bwd_s_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 5e-3, "element " << t << " must be approx same");
    }
}

void expl_black_scholes_equation_dirichlet_bc_euler()
{
    using lss_pde_solvers::default_heat_solver_configs::host_expl_bwd_euler_solver_config_ptr;
    using lss_pde_solvers::one_dimensional::explicit_solvers::heat_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Black-Scholes Call equation: \n\n";
    std::cout << " Using Euler method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(x,t) = 0.5*sig*sig*x*x*U_xx(x,t) + r*x*U_x(x,t) - "
                 "r*U(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < x < 20 and 0 < t < 1,\n";
    std::cout << " U(0,t) = 0 and  U(20,t) = 20-K*exp(-r*(1-t)),0 < t < 1 \n\n";
    std::cout << " U(x,T) = max(0,x-K), x in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10;
    auto const &maturity = 1.0;
    auto const &rate = 0.2;
    auto const &sig = 0.25;
    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 6000;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 20.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_1d>(space_range, Sd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double x) { return 0.5 * sig * sig * x * x; };
    auto b = [=](double t, double x) { return rate * x; };
    auto c = [=](double t, double x) { return -rate; };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_1d>(a, b, c);
    // terminal condition:
    auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_bwd_euler_solver_config_ptr);
    // prepare container for solution:
    container_t solution(double{}, Sd);
    // get the solution:
    pdesolver.solve(solution);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    std::cout << "tp : FDM | Exact | Abs Diff\n";
    double benchmark{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark = bs_exact.call(x);
        std::cout << "t_" << j << ": " << solution[j] << " |  " << benchmark << " | " << (solution[j] - benchmark)
                  << '\n';
    }
    // TEST:
    for (std::size_t t = 0; t < solution.size(); ++t)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, t));
        benchmark = bs_exact.call(x);
        LSS_ASSERT(std::abs(solution[t] - benchmark) < 3e-3, "element " << t << " must be approx same");
    }
}

void test_expl_black_scholes_equation_dirichlet_bc_ade()
{
    std::cout << "============================================================\n";
    std::cout << "========= Explicit Black-Scholes Equation (Dir BC) =========\n";
    std::cout << "============================================================\n";

    expl_black_scholes_equation_dirichlet_bc_barakat_clark();
    expl_black_scholes_equation_dirichlet_bc_saulyev();
    expl_black_scholes_equation_dirichlet_bc_euler();

    std::cout << "============================================================\n";
}
