#include "lss_xml_t.hpp"

#include <cmath>
#include <sstream>

#include "../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../boundaries/lss_neumann_boundary.hpp"
#include "../../boundaries/lss_robin_boundary.hpp"
#include "../../common/lss_xml.hpp"
#include "../../containers/lss_matrix_2d.hpp"
#include "../../containers/lss_matrix_3d.hpp"
#include "../../ode_solvers/second_degree/lss_ode_equation.hpp"
#include "../../pde_solvers/1d/heat_type/lss_heat_equation.hpp"
#include "../../pde_solvers/1d/wave_type/lss_wave_equation.hpp"
#include "../../pde_solvers/2d/heat_type/lss_heston_equation.hpp"

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::dirichlet_boundary_2d;
using lss_boundary::dirichlet_boundary_3d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::neumann_boundary_2d;
using lss_boundary::neumann_boundary_3d;
using lss_boundary::robin_boundary_1d;
using lss_containers::matrix_2d;
using lss_containers::matrix_3d;
using lss_enumerations::grid_enum;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_enumerations::splitting_method_enum;
using lss_grids::grid_1d;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d;
using lss_grids::grid_config_hints_2d;
using lss_grids::grid_config_hints_3d;
using lss_grids::grid_transform_config_1d;
using lss_ode_solvers::ode_coefficient_data_config;
using lss_ode_solvers::ode_data_config;
using lss_ode_solvers::ode_discretization_config;
using lss_ode_solvers::ode_nonhom_data_config;
using lss_ode_solvers::implicit_solvers::ode_equation;
using lss_pde_solvers::heat_coefficient_data_config_1d;
using lss_pde_solvers::heat_coefficient_data_config_2d;
using lss_pde_solvers::heat_coefficient_data_config_3d;
using lss_pde_solvers::heat_data_config_1d;
using lss_pde_solvers::heat_data_config_2d;
using lss_pde_solvers::heat_data_config_3d;
using lss_pde_solvers::heat_implicit_solver_config;
using lss_pde_solvers::heat_initial_data_config_1d;
using lss_pde_solvers::heat_initial_data_config_2d;
using lss_pde_solvers::heat_initial_data_config_3d;
using lss_pde_solvers::pde_discretization_config_1d;
using lss_pde_solvers::pde_discretization_config_2d;
using lss_pde_solvers::pde_discretization_config_3d;
using lss_pde_solvers::splitting_method_config;
using lss_pde_solvers::wave_coefficient_data_config_1d;
using lss_pde_solvers::wave_data_config_1d;
using lss_pde_solvers::wave_implicit_solver_config;
using lss_pde_solvers::wave_initial_data_config_1d;
using lss_pde_solvers::one_dimensional::implicit_solvers::heat_equation;
using lss_pde_solvers::one_dimensional::implicit_solvers::wave_equation;
using lss_pde_solvers::two_dimensional::implicit_solvers::heston_equation;
using lss_utility::black_scholes_exact;
using lss_utility::container_t;
using lss_utility::pi;
using lss_utility::range;
using lss_xml::xml;

// ODEs

void impl_simple_ode_thomes_lu_qr_xml()
{
    using lss_ode_solvers::default_ode_solver_configs::dev_cusolver_qr_solver_config_ptr;
    const std::string test_name{"ode_bvp_neumann_robin"};

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
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
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
    auto const &neumann = [](double t) { return -1.0; };
    auto const &robin_first = [](double t) { return 2.0; };
    auto const &robin_second = [](double t) { return 0.0; };
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

    // print both of these
    std::stringstream ssa;
    ssa << "outputs/xmls/" << test_name << ".xml";
    std::string file_name_approx{ssa.str()};
    std::ofstream approx(file_name_approx);
    xml(discretization_ptr, grid_config_hints_ptr, solution, approx);
    approx.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_approx << "\n";
}

void test_impl_simple_ode_thomes_lu_qr_xml()
{
    std::cout << "============================================================\n";
    std::cout << "=========== Implicit Simple ODE (Thomas LU)  ===============\n";
    std::cout << "============================================================\n";

    impl_simple_ode_thomes_lu_qr_xml();

    std::cout << "============================================================\n";
}

// PDEs
void impl_bs_thomas_lu_solver_euler_crv_xml()
{
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_euler_solver_config_ptr;
    const std::string test_name_bench{"bs_bvp_thomas_lu_solver_euler_bench"};
    const std::string test_name_num{"bs_bvp_thomas_lu_solver_euler_numerical"};

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
    auto terminal_condition = [=](double x) { return std::max(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid hints:
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
    // get the benchmark:
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    container_t benchmark(solution.size());
    double x{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark[j] = bs_exact.call(x);
    }
    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solution, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_bs_thomas_lu_solver_cn_crv_xml()
{
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;
    const std::string test_name_bench{"bs_bvp_thomas_lu_solver_cn_bench"};
    const std::string test_name_num{"bs_bvp_thomas_lu_solver_cn_numerical"};

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
    auto terminal_condition = [=](double x) { return std::max(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid hints:
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

    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    container_t benchmark(solution.size());
    double x{};
    for (std::size_t j = 0; j < solution.size(); ++j)
    {
        x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
        benchmark[j] = bs_exact.call(x);
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solution, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_bs_thomas_lu_crv_xml()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Black-Scholes (Thomas LU) Equation (Dir BC) ====\n";
    std::cout << "============================================================\n";

    impl_bs_thomas_lu_solver_euler_crv_xml();
    impl_bs_thomas_lu_solver_cn_crv_xml();

    std::cout << "============================================================\n";
}

void impl_bs_thomas_lu_euler_srf_xml()
{
    const std::string test_name_bench{"bs_bvp_thomas_lu_solver_euler_surf_bench"};
    const std::string test_name_num{"bs_bvp_thomas_lu_solver_euler_surf_numerical"};

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
    auto terminal_condition = [=](double x) { return std::max(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid hints:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(strike);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_tlusolver_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    // get the benchmark:
    double x{};
    double const k = discretization_ptr->time_step();
    matrix_2d benchmark(Td, Sd);
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, bs_exact.call(x, maturity - t * k));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_bs_thomas_lu_cn_srf_xml()
{
    const std::string test_name_bench{"bs_bvp_thomas_lu_solver_cn_surf_bench"};
    const std::string test_name_num{"bs_bvp_thomas_lu_solver_cn_surf_numerical"};
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
    auto terminal_condition = [=](double x) { return std::max(0.0, x - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_1d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_1d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // boundary conditions:
    auto const &dirichlet_low = [=](double t) { return 0.0; };
    auto const &dirichlet_high = [=](double t) { return (20.0 - strike * std::exp(-rate * (maturity - t))); };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_1d>(dirichlet_high);
    auto const &boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_1d>(strike, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    black_scholes_exact bs_exact(0.0, strike, rate, sig, maturity);

    double x{};
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, bs_exact.call(x, maturity - t * k));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_bs_thomas_lu_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "== Implicit Black-Scholes (Thomas LU) Equation (Dir BC) ====\n";
    std::cout << "============================================================\n";

    impl_bs_thomas_lu_euler_srf_xml();
    impl_bs_thomas_lu_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_ph_dirichlet_bvp_device_qr_euler_srf_xml()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_euler_solver_config_ptr;
    const std::string test_name_bench{"pure_heat_bvp_cuda_solver_device_qr_euler_surf_bench"};
    const std::string test_name_num{"pure_heat_bvp_cuda_solver_device_qr_euler_surf_numerical"};

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
    // grid hints:
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
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k, 20));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_ph_dirichlet_bvp_device_qr_cn_srf_xml()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_fwd_cusolver_qr_cn_solver_config_ptr;
    const std::string test_name_bench{"pure_heat_bvp_cuda_solver_device_qr_cn_surf_bench"};
    const std::string test_name_num{"pure_heat_bvp_cuda_solver_device_qr_cn_surf_numerical"};

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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
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
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k, 20));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_ph_dirichlet_bvp_device_qr_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << " Implicit Pure Heat (CUDA QR DEVICE) Equation (Dirichlet BC)\n";
    std::cout << "============================================================\n";

    impl_ph_dirichlet_bvp_device_qr_euler_srf_xml();
    impl_ph_dirichlet_bvp_device_qr_cn_srf_xml();

    std::cout << "============================================================\n";
}

void expl_ph_neumann_neumann_bvp_euler_srf_xml()
{
    const std::string test_name_bench{"pure_heat_bvp_neumann_euler_surf_bench"};
    const std::string test_name_num{"pure_heat_bvp_neumann_euler_surf_numerical"};

    using lss_pde_solvers::default_heat_solver_configs::dev_expl_fwd_euler_solver_config_ptr;
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
    std::size_t const Sd = 50;
    // number of time subdivisions:
    std::size_t const Td = 2000;
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
    // grid hints:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>();
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_expl_fwd_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        double const pipi = pi() * pi();
        double const first = 4.0 / pipi;
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
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k, 20));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_expl_ph_neumann_bvp_euler_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "===== Explicit Pure Heat (ADE ) Equation (Non-Dir BC) ======\n";
    std::cout << "============================================================\n";

    expl_ph_neumann_neumann_bvp_euler_srf_xml();

    std::cout << "============================================================\n";
}

void impl_adv_thomas_lu_solver_cn_srf_xml()
{
    const std::string test_name_bench{"adv_bvp_thomas_lu_surf_bench"};
    const std::string test_name_num{"adv_bvp_thomas_lu_surf_numerical"};
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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    heat_equation pdesolver(heat_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
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
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k, 20));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_adv_thomas_lu_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "========= Implicit Advection Equation (Non-Dir BC) =========\n";
    std::cout << "============================================================\n";

    impl_adv_thomas_lu_solver_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_pw_dirichlet_bvp_cuda_device_qr_srf_xml()
{
    const std::string test_name_bench{"adv_bvp_thomas_lu_surf_bench"};
    const std::string test_name_num{"adv_bvp_thomas_lu_surf_numerical"};
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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t, std::size_t n) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
    };

    double x{};
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k, 20));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_pw_dirichlet_bvp_cuda_device_qr_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "========= Implicit Pure Wave Equation (Non-Dir BC) =========\n";
    std::cout << "============================================================\n";

    impl_pw_dirichlet_bvp_cuda_device_qr_srf_xml();

    std::cout << "============================================================\n";
}

void impl_w_dirichlet_bvp_host_lu_srf_xml()
{
    const std::string test_name_bench{"w_bvp_thomas_lu_srf_bench"};
    const std::string test_name_num{"w_bvp_thomas_lu_srf_numerical"};
    using lss_pde_solvers::default_wave_solver_configs::host_fwd_tlusolver_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Wave equation: \n\n";
    std::cout << " Using Thomas LU Solver method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_tt(x,t) = 4*U_xx(x,t), \n\n";
    std::cout << " where\n\n";
    std::cout << " x in <0,1> and t > 0,\n";
    std::cout << " U(0,t) = -4*t*t,U(1,t) = -4*t*t + 8*t, t > 0 \n\n";
    std::cout << " U(x,0) = x*(1 - x), x in <0,1> \n\n";
    std::cout << " U_x(x,0) = 0, x in <0,1> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 5.7);
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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_tlusolver_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double res = x - x * x - 4.0 * t * t + 8.0 * t * x;
        return (res);
    };

    double x{};
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_w_dirichlet_bvp_host_lu_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "===== Implicit Wave (Thomas LU) Equation (Dirichlet BC) ====\n";
    std::cout << "============================================================\n";

    impl_w_dirichlet_bvp_host_lu_srf_xml();

    std::cout << "============================================================\n";
}

void impl_w_bvp_host_double_sweep_srf_xml()
{
    const std::string test_name_bench{"w_bvp_thomas_dss_srf_bench"};
    const std::string test_name_num{"w_bvp_thomas_dss_srf_numerical"};
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
    auto const &space_range = std::make_shared<range>(0.0, pi());
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 5.7);
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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_fwd_dssolver_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double exp_half = std::exp(-0.5 * t);
        const double sqrt_3 = std::sqrt(3.0);
        const double arg = 0.5 * sqrt_3;
        const double res = exp_half * std::sin(x) * (std::cos(arg * t) + (sin(arg * t) / sqrt_3));
        return (res);
    };

    double x{};
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_w_bvp_host_dss_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "=== Implicit Wave (Double Sweep) Equation (Dirichlet BC) ===\n";
    std::cout << "============================================================\n";

    impl_w_bvp_host_double_sweep_srf_xml();

    std::cout << "============================================================\n";
}

void impl_pw_neumann_bvp_cuda_device_qr_srf_xml()
{
    const std::string test_name_bench{"pw_bvp_cuda_device_qr_srf_bench"};
    const std::string test_name_num{"pw_bvp_cuda_device_qr_srf_numerical"};
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
    std::cout << " U_x(x,0) = 1 - cos(4*x), x in <0,pi> \n\n";
    std::cout << "============================================================\n";

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range = std::make_shared<range>(0.0, pi());
    // time range
    auto const &time_range = std::make_shared<range>(0.0, 3.8);
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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            dev_fwd_cusolver_qr_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double var1 = 3.0 * std::cos(2.0 * t) * std::cos(x);
        const double var2 = -0.125 * std::sin(8.0 * t) * std::cos(4.0 * x);
        return (t + var1 + var2);
    };

    double x{};
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_pw_neumann_bvp_cuda_device_qr_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Wave (CUDA QR DEVICE) Equation (Neumann BC)=\n";
    std::cout << "============================================================\n";

    impl_pw_neumann_bvp_cuda_device_qr_srf_xml();

    std::cout << "============================================================\n";
}

void expl_pw_dirichlet_bvp_cuda_host_srf_xml()
{
    const std::string test_name_bench{"pw_expl_bvp_cuda_host_qr_srf_bench"};
    const std::string test_name_num{"pw_expl_bvp_cuda_host_qr_srf_numerical"};
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
    // grid hints:
    auto const alpha_scale = 3.;
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_1d>(0.5, alpha_scale, grid_enum::Nonuniform);
    // initialize pde solver
    wave_equation pdesolver(wave_data_ptr, discretization_ptr, boundary_pair, grid_config_hints_ptr,
                            host_expl_fwd_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solutions(Td, Sd);
    // get the solution:
    pdesolver.solve(solutions);
    // get exact solution:
    auto exact = [](double x, double t) {
        const double var1 = std::sin(pi() * x);
        const double var2 = std::cos(pi() * t);
        return (var1 * var2);
    };

    double x{};
    double const k = discretization_ptr->time_step();
    auto const grid_cfg = std::make_shared<grid_config_1d>(discretization_ptr);
    auto const grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_ptr, grid_config_hints_ptr);
    matrix_2d benchmark(Td, Sd);
    for (std::size_t t = 0; t < solutions.rows(); ++t)
    {
        for (std::size_t j = 0; j < solutions.columns(); ++j)
        {
            x = grid_1d::transformed_value(grid_trans_cfg, grid_1d::value(grid_cfg, j));
            benchmark(t, j, exact(x, t * k));
        }
    }

    // xml both of these
    std::stringstream ssb;
    ssb << "outputs/xmls/" << test_name_bench << ".xml";
    std::string file_name_bench{ssb.str()};
    std::ofstream bench(file_name_bench);
    xml(discretization_ptr, grid_config_hints_ptr, solutions, bench);
    bench.close();
    std::cout << "Explicit solution has been saved to XML file: " << file_name_bench << "\n";
    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, benchmark, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_expl_pw_dirichlet_bvp_cuda_host_host_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "= Implicit Pure Wave (CUDA DEVICE) Equation (Dirichlet BC) =\n";
    std::cout << "============================================================\n";

    expl_pw_dirichlet_bvp_cuda_host_srf_xml();

    std::cout << "============================================================\n";
}

void impl_heston_cuda_qr_crank_nicolson_dr_srf_xml()
{
    const std::string test_name_num{"heston_cuda_qr_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using CUDA QR algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < s < 20, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(0,v,t) = 0 and  U_s(20,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 10.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 100;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 150;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(0.0, 20.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha = 3.;
    auto const beta = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha, beta, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, dev_bwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_heston_cuda_qr_cn_dr_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "=========== Implicit Heston Equation (CUDA DEVICE) =========\n";
    std::cout << "============================================================\n";

    impl_heston_cuda_qr_crank_nicolson_dr_srf_xml();

    std::cout << "============================================================\n";
}

void impl_heston_thomas_lu_cn_srf_xml()
{
    const std::string test_name_num{"heston_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(30.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_heston_upoutcall_barrier_thomas_lu_cn_srf_xml()
{
    const std::string test_name_num{"heston_upoutcall_barrier_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " B = 150, 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U(B,v,t) = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U_v(s,1,t) = 0, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K) if s < B, else 0.0 for s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &barrier = 150.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return ((s < barrier) ? std::max(0.0, s - strike) : 0.0); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_s_low = [=](double t, double v) { return 0.0; };
    auto const &dirichlet_s_high = [=](double t, double v) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &neumann_high = [=](double t, double s) { return 0.0; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_heston_downoutcall_barrier_put_thomas_lu_cn_srf_xml()
{
    const std::string test_name_num{"heston_downoutcall_barrier_put_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " B = 150, 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U(B,v,t) = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U_v(s,1,t) = 0, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,K - s) if s < B, else 0.0 for s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &barrier = 80.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return ((s > barrier) ? std::max(0.0, s - strike) : 0.0); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_s_low = [=](double t, double v) { return 0.0; };
    auto const &dirichlet_s_high = [=](double t, double v) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &neumann_high = [=](double t, double s) { return 0.0; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_heston_upoutput_barrier_put_thomas_lu_cn_srf_xml()
{
    const std::string test_name_num{"heston_upoutput_barrier_put_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " B = 150, 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U(B,v,t) = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U_v(s,1,t) = 0, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,K - s) if s < B, else 0.0 for s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &barrier = 150.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return ((s < barrier) ? std::max(0.0, strike - s) : 0.0); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_s_low = [=](double t, double v) { return 0.0; };
    auto const &dirichlet_s_high = [=](double t, double v) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &neumann_high = [=](double t, double s) { return 0.0; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_heston_downoutput_barrier_put_thomas_lu_cn_srf_xml()
{
    const std::string test_name_num{"heston_downoutput_barrier_put_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " B = 150, 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U(B,v,t) = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U_v(s,1,t) = 0, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,K - s) if s < B, else 0.0 for s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &barrier = 80.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return ((s > barrier) ? std::max(0.0, strike - s) : 0.0); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_s_low = [=](double t, double v) { return 0.0; };
    auto const &dirichlet_s_high = [=](double t, double v) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &neumann_high = [=](double t, double s) { return 0.0; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_heston_thomas_lu_cn_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "======== Implicit Heston Equation (Thomas LU Solver) =======\n";
    std::cout << "============================================================\n";

    impl_heston_thomas_lu_cn_srf_xml();
    impl_heston_upoutcall_barrier_thomas_lu_cn_srf_xml();
    impl_heston_downoutcall_barrier_put_thomas_lu_cn_srf_xml();
    impl_heston_upoutput_barrier_put_thomas_lu_cn_srf_xml();
    impl_heston_downoutput_barrier_put_thomas_lu_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_sabr_double_sweep_cn_srf_xml()
{
    const std::string test_name_num{"sabr_double_lu_cn_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value SABR Call equation: \n\n";
    std::cout << " Using Double Sweep algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,a,t) = 0.5*a*a*s^(2b)*D^(2*(1-b))*U_ss(s,a,t) "
                 " + 0.5*sig*sig*a*a*U_vv(s,v,t)"
                 " + rho*sig*s^b*D^(1-b)*a*a*U_sv(s,a,t) + r*s*U_s(s,a,t)"
                 " - r*U(s,a,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(0,a,t) = 0 and  U_s(200,a,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t) - rU(s,0,t) - U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,a,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.081;
    auto const &rho = 0.6;
    auto const &beta = 0.7;
    // number of space subdivisions for spot:
    std::size_t const Sd = 100;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto D = [=](double t, double s, double alpha) { return std::exp(-rate * (maturity - t)); };
    auto a = [=](double t, double s, double alpha) {
        return (0.5 * alpha * alpha * std::pow(s, 2.0 * beta) * std::pow(D(t, s, alpha), 2.0 * (1.0 - beta)));
    };
    auto b = [=](double t, double s, double alpha) { return (0.5 * sig_sig * sig_sig * alpha * alpha); };
    auto c = [=](double t, double s, double alpha) {
        return (rho * sig_sig * alpha * alpha * std::pow(s, beta) * std::pow(D(t, s, alpha), (1.0 - beta)));
    };
    auto d = [=](double t, double s, double alpha) { return (rate * s); };
    auto e = [=](double t, double s, double alpha) { return 0.0; };
    auto f = [=](double t, double s, double alpha) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_sabr_barrier_double_sweep_cn_srf_xml()
{
    const std::string test_name_num{"sabr_barrier_double_lu_cn_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value SABR Barrier up-and-call equation: \n\n";
    std::cout << " Using Double Sweep algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,a,t) = 0.5*a*a*s^(2b)*D^(2*(1-b))*U_ss(s,a,t) "
                 " + 0.5*sig*sig*a*a*U_vv(s,v,t)"
                 " + rho*sig*s^b*D^(1-b)*a*a*U_sv(s,a,t) + r*s*U_s(s,a,t)"
                 " - r*U(s,a,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " B = 150, 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(0,a,t) = 0 and  U(B,a,t) = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t) - rU(s,0,t) - U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U_v(s,1,t) = 0, 0 < t < 1\n";
    std::cout << " U(s,a,T) = max(0,s - K) s < B else 0, s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &barrier = 150.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.081;
    auto const &rho = 0.6;
    auto const &beta = 0.7;
    // number of space subdivisions for spot:
    std::size_t const Sd = 100;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto D = [=](double t, double s, double alpha) { return std::exp(-rate * (maturity - t)); };
    auto a = [=](double t, double s, double alpha) {
        return (0.5 * alpha * alpha * std::pow(s, 2.0 * beta) * std::pow(D(t, s, alpha), 2.0 * (1.0 - beta)));
    };
    auto b = [=](double t, double s, double alpha) { return (0.5 * sig_sig * sig_sig * alpha * alpha); };
    auto c = [=](double t, double s, double alpha) {
        return (rho * sig_sig * alpha * alpha * std::pow(s, beta) * std::pow(D(t, s, alpha), (1.0 - beta)));
    };
    auto d = [=](double t, double s, double alpha) { return (rate * s); };
    auto e = [=](double t, double s, double alpha) { return 0.0; };
    auto f = [=](double t, double s, double alpha) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return ((s < barrier) ? std::max(0.0, s - strike) : 0.0); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_s_low = [=](double t, double v) { return 0.0; };
    auto const &dirichlet_s_high = [=](double t, double v) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &neumann_high = [=](double t, double s) { return 0.0; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_sabr_double_sweep_cn_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "====== Implicit SABR Equation (Double Sweep Solver) ========\n";
    std::cout << "============================================================\n";

    impl_sabr_double_sweep_cn_srf_xml();
    impl_sabr_barrier_double_sweep_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_heston_thomas_lu_dr_cn_srf_xml()
{
    const std::string test_name_num{"heston_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 150.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_2d>();

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_heston_thomas_lu_dr_cn_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "===== Implicit Heston Equation (Thomas LU Solver => DR) ====\n";
    std::cout << "============================================================\n";

    impl_heston_thomas_lu_dr_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_heston_dss_solver_cs_cn_srf_xml()
{
    const std::string test_name_num{"heston_thomas_lu_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr = std::make_shared<splitting_method_config>(splitting_method_enum::CraigSneyd);
    // grid:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void impl_heston_put_dss_solver_cs_cn_srf_xml()
{
    const std::string test_name_num{"heston_put_dss_cn_dr_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Put equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = K*exp{-r(T-t)} and  U(200,v,t) = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = K*exp{-r(T-t)}, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,K - s), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, strike - s); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_s_low = [=](double t, double v) { return strike * std::exp(-rate * (maturity - t)); };
    auto const &dirichlet_s_high = [=](double t, double v) { return 0.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_low);
    auto const &boundary_high_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_s_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return strike * std::exp(-rate * (maturity - t)); };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr = std::make_shared<splitting_method_config>(splitting_method_enum::CraigSneyd);
    // grid:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_heston_dss_cs_cn_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "===== Implicit Heston Equation (Thomas LU Solver => CS) ====\n";
    std::cout << "============================================================\n";

    // impl_heston_dss_solver_cs_cn_srf_xml();
    impl_heston_put_dss_solver_cs_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_heston_thomas_lu_mcs_cn_srf_xml()
{
    const std::string test_name_num{"heston_thomas_lu_mcs_cn_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 180.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::ModifiedCraigSneyd);
    // grid:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_heston_thomas_lu_mcs_cn_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "==== Implicit Heston Equation (Thomas LU Solver => MCS) ====\n";
    std::cout << "============================================================\n";

    impl_heston_thomas_lu_mcs_cn_srf_xml();

    std::cout << "============================================================\n";
}

void impl_heston_thomas_lu_hw_cn_srf_xml()
{
    const std::string test_name_num{"heston_thomas_lu_hw_cn_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.2;
    // number of space subdivisions for spot:
    std::size_t const Sd = 150;
    // number of space subdivision for volatility:
    std::size_t const Vd = 50;
    // number of time subdivisions:
    std::size_t const Td = 200;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::HundsdorferVerwer);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_impl_heston_thomas_lu_hw_cn_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "===== Implicit Heston Equation (Thomas LU Solver => HV) ====\n";
    std::cout << "============================================================\n";

    impl_heston_thomas_lu_hw_cn_srf_xml();

    std::cout << "============================================================\n";
}

// void impl_hhw_dsssolver_dr_cn_srf_xml()
//{
//     const std::string test_name_num{"hhw_dsssolver_dr_cn_srf_numerical"};
//     using lss_pde_solvers::default_heat_solver_configs::host_bwd_dssolver_cn_solver_config_ptr;
//
//     std::cout << "============================================================\n";
//     std::cout << "Solving Boundary-value Heston-Hull-White Call equation: \n\n";
//     std::cout << "Using DSS algo with implicit Crank-Nicolson method\n\n";
//     std::cout << " Value type: " << typeid(double).name() << "\n\n";
//     std::cout << " U_t(t,s,v,r) = 0.5*s*s*v*U_ss(t,s,v,r) +"
//                  "0.5*vol_1*vol_1*v*U_vv(t,s,v,r) +"
//                  "0.5*vol_2*vol_2*U_rr(t,s,v,r) + "
//                  "rho_12*vol_1*s*v*U_sv(t,s,v,r) + "
//                  "rho_13*vol_2*s*sqrt(v)*U_sr(t,s,v,r) + "
//                  "rho_23*vol_1*vol_2*sqrt(v)*U_vr(t,s,v,r) + "
//                  "r*s*U_s(t,s,v,r) + [k*(theta-v)-lambda*v]*U_v(t,s,v,r) +"
//                  "a*(b-r)*U_r(t,s,v,r) - r*U(t,s,v,r)\n\n";
//     std::cout << " where\n\n";
//     std::cout << " 0 < s < 20, 0 < v < 1,-1 < r < 1, and 0 < t < 1,\n";
//     std::cout << " U(t,0,v,r) = 0 and  U_s(t,20,v,r) - 1 = 0, 0 < t < 1\n";
//     std::cout << " 0.5*vol_2*vol_2*U_rr(t,s,0,r) + r*s*U_s(t,s,0,r) + "
//                  "k*theta*U_v(t,s,0,r) + a*(b-r)*U_r(t,s,0,r) - "
//                  "rU(t,s,0,r) - U_t(t,s,0,r) = 0,"
//                  "0 < t < 1\n";
//     std::cout << " U(t,s,1,r) = s, 0 < t < 1\n";
//     std::cout << " U_r(t,s,v,-1) = 0 and  U_r(t,s,v,1) = 0, 0 < t < 1\n";
//     std::cout << " U(T,s,v,r) = max(0,s - K), s in <0,20> \n\n";
//     std::cout << "============================================================\n";
//
//     // set up call option parameters:
//     auto const &strike = 100.0;
//     auto const &maturity = 1.0;
//     auto const &v_sig = 0.7;
//     auto const &v_kappa = 3.0;
//     auto const &v_theta = 0.2;
//     auto const &rho_12 = 0.6;
//     auto const &rho_13 = 0.2;
//     auto const &rho_23 = 0.4;
//     auto const &c_1 = 0.3;
//     auto const &c_2 = 0.01;
//     auto const &c_3 = 0.02;
//     auto const &r_a = 0.2;
//     auto const &r_sig = 0.1;
//     // number of space subdivisions for spot:
//     std::size_t const Sd = 55;
//     // number of space subdivision for volatility:
//     std::size_t const Vd = 28;
//     // number of space subdivision for rate:
//     std::size_t const Rd = 28;
//     // number of time subdivisions:
//     std::size_t const Td = 500;
//     // space Spot range:
//     auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
//     // space Vol range:
//     auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
//     // space Rate range:
//     auto const &spacez_range = std::make_shared<range>(0.12, 0.21);
//     // time range
//     auto const &time_range = std::make_shared<range>(0.0, maturity);
//     // discretization config:
//     auto const discretization_ptr = std::make_shared<pde_discretization_config_3d>(
//         spacex_range, spacey_range, spacez_range, Sd, Vd, Rd, time_range, Td);
//     // coeffs:
//     auto B = [=](double t) { return (c_1 - c_2 * std::exp(-c_3 * (maturity - t))); };
//     auto a = [=](double t, double s, double v, double r) { return (0.5 * v * s * s); };
//     auto b = [=](double t, double s, double v, double r) { return (0.5 * v_sig * v_sig * v); };
//     auto c = [=](double t, double s, double v, double r) { return (0.5 * r_sig * r_sig); };
//     auto d = [=](double t, double s, double v, double r) { return (rho_12 * v_sig * s * v); };
//     auto e = [=](double t, double s, double v, double r) { return (rho_13 * r_sig * s * sqrt(v)); };
//     auto f = [=](double t, double s, double v, double r) { return (rho_23 * v_sig * r_sig * sqrt(v)); };
//     auto g = [=](double t, double s, double v, double r) { return (r * s); };
//     auto h = [=](double t, double s, double v, double r) { return (v_kappa * (v_theta - v)); };
//     auto i = [=](double t, double s, double v, double r) { return (r_a * (B(t) - r)); };
//     auto j = [=](double t, double s, double v, double r) { return (-r); };
//     auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_3d>(a, b, c, d, e, f, g, h, i,
//     j);
//     // terminal condition:
//     auto terminal_condition = [=](double s, double v, double r) { return std::max(0.0, s - strike); };
//     auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_3d>(terminal_condition);
//     // heat data config:
//     auto const heat_data_ptr = std::make_shared<heat_data_config_3d>(heat_coeffs_data_ptr, heat_init_data_ptr);
//     // spot boundary conditions:
//     auto const &dirichlet_s = [=](double t, double v, double r) { return 0.0; };
//     auto const &neumann_s = [=](double t, double v, double r) { return -1.0; };
//     auto const &boundary_low_s_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_s);
//     auto const &boundary_high_s_ptr = std::make_shared<neumann_boundary_3d>(neumann_s);
//     auto const &x_boundary_pair = std::make_pair(boundary_low_s_ptr, boundary_high_s_ptr);
//     //  upper vol boundary:
//     auto const &dirichlet_v = [=](double t, double s, double r) { return s; };
//     auto const &y_upper_boundary_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_v);
//     // rate boundary conditions:
//     auto const &neumann_low_r = [=](double t, double s, double v) { return 0.0; };
//     auto const &neumann_high_r = [=](double t, double s, double v) { return 0.0; };
//     auto const &boundary_low_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_low_r);
//     auto const &boundary_high_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_high_r);
//     auto const &z_boundary_pair = std::make_pair(boundary_low_r_ptr, boundary_high_r_ptr);
//     // splitting method configuration:
//     auto const &splitting_config_ptr =
//         std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
//
//     // grid config:
//     auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_3d>(strike);
//
//     // initialize pde solver
//     hhw_equation pdesolver(heat_data_ptr, discretization_ptr, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair,
//                            splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_cn_solver_config_ptr);
//     // prepare container for solution:
//     container_3d<by_enum::RowPlane> solution(Sd, Vd, Rd, double{});
//     // get the solution:
//     pdesolver.solve(solution);
//
//     std::stringstream ssn;
//     ssn << "outputs/xmls/" << test_name_num << ".xml";
//     std::string file_name_num{ssn.str()};
//     std::ofstream numerical(file_name_num);
//     xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
//     numerical.close();
//     std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
// }

// void impl_hhw_dsssolver_dr_0_66_cn_srf_xml()
//{
//     const std::string test_name_num{"hhw_dsssolver_dr_0_66_cn_srf_numerical"};
//     using lss_pde_solvers::default_heat_solver_configs::build_implicit_config;
//
//     auto const host_bwd_dssolver_0_66_cn_solver_config_ptr =
//         build_implicit_config(memory_space_enum::Host, traverse_direction_enum::Backward,
//                               tridiagonal_method_enum::ThomasLUSolver, factorization_enum::None, (2.0 / 3.0));
//
//     std::cout << "============================================================\n";
//     std::cout << "Solving Boundary-value Heston-Hull-White Call equation: \n\n";
//     std::cout << "Using DSS algo with implicit Crank-Nicolson method\n\n";
//     std::cout << " Value type: " << typeid(double).name() << "\n\n";
//     std::cout << " U_t(t,s,v,r) = 0.5*s*s*v*U_ss(t,s,v,r) +"
//                  "0.5*vol_1*vol_1*v*U_vv(t,s,v,r) +"
//                  "0.5*vol_2*vol_2*U_rr(t,s,v,r) + "
//                  "rho_12*vol_1*s*v*U_sv(t,s,v,r) + "
//                  "rho_13*vol_2*s*sqrt(v)*U_sr(t,s,v,r) + "
//                  "rho_23*vol_1*vol_2*sqrt(v)*U_vr(t,s,v,r) + "
//                  "r*s*U_s(t,s,v,r) + [k*(theta-v)-lambda*v]*U_v(t,s,v,r) +"
//                  "a*(b-r)*U_r(t,s,v,r) - r*U(t,s,v,r)\n\n";
//     std::cout << " where\n\n";
//     std::cout << " 0 < s < 20, 0 < v < 1,-1 < r < 1, and 0 < t < 1,\n";
//     std::cout << " U(t,0,v,r) = 0 and  U_s(t,20,v,r) - 1 = 0, 0 < t < 1\n";
//     std::cout << " 0.5*vol_2*vol_2*U_rr(t,s,0,r) + r*s*U_s(t,s,0,r) + "
//                  "k*theta*U_v(t,s,0,r) + a*(b-r)*U_r(t,s,0,r) - "
//                  "rU(t,s,0,r) - U_t(t,s,0,r) = 0,"
//                  "0 < t < 1\n";
//     std::cout << " U(t,s,1,r) = s, 0 < t < 1\n";
//     std::cout << " U_r(t,s,v,-1) = 0 and  U_r(t,s,v,1) = 0, 0 < t < 1\n";
//     std::cout << " U(T,s,v,r) = max(0,s - K), s in <0,20> \n\n";
//     std::cout << "============================================================\n";
//
//     // set up call option parameters:
//     auto const &strike = 100.0;
//     auto const &maturity = 15;
//     auto const &v_sig = 0.9;
//     auto const &v_kappa = 0.3;
//     auto const &v_theta = 0.04;
//     auto const &rho_12 = -0.5;
//     auto const &rho_13 = 0.2;
//     auto const &rho_23 = 0.1;
//     auto const &c_1 = 0.055;
//     auto const &c_2 = 0.025;
//     auto const &c_3 = 1.6;
//     auto const &r_a = 0.16;
//     auto const &r_sig = 0.3;
//
//     // number of space subdivisions for spot:
//     std::size_t const Sd = 50;
//     // number of space subdivision for volatility:
//     std::size_t const Vd = 25;
//     // number of space subdivision for rate:
//     std::size_t const Rd = 24;
//     // number of time subdivisions:
//     std::size_t const Td = 500;
//     // space Spot range:
//     auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
//     // space Vol range:
//     auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
//     // space Rate range:
//     auto const &spacez_range = std::make_shared<range>(0.12, 0.5);
//     // time range
//     auto const &time_range = std::make_shared<range>(0.0, maturity);
//     // discretization config:
//     auto const discretization_ptr = std::make_shared<pde_discretization_config_3d>(
//         spacex_range, spacey_range, spacez_range, Sd, Vd, Rd, time_range, Td);
//     // coeffs:
//     auto B = [=](double t) { return (c_1 - c_2 * std::exp(-c_3 * (maturity - t))); };
//     auto a = [=](double t, double s, double v, double r) { return (0.5 * v * s * s); };
//     auto b = [=](double t, double s, double v, double r) { return (0.5 * v_sig * v_sig * v); };
//     auto c = [=](double t, double s, double v, double r) { return (0.5 * r_sig * r_sig); };
//     auto d = [=](double t, double s, double v, double r) { return (rho_12 * v_sig * s * v); };
//     auto e = [=](double t, double s, double v, double r) { return (rho_13 * r_sig * s * sqrt(v)); };
//     auto f = [=](double t, double s, double v, double r) { return (rho_23 * v_sig * r_sig * sqrt(v)); };
//     auto g = [=](double t, double s, double v, double r) { return (r * s); };
//     auto h = [=](double t, double s, double v, double r) { return (v_kappa * (v_theta - v)); };
//     auto i = [=](double t, double s, double v, double r) { return (r_a * (B(t) - r)); };
//     auto j = [=](double t, double s, double v, double r) { return (-r); };
//     auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_3d>(a, b, c, d, e, f, g, h, i,
//     j);
//     // terminal condition:
//     auto terminal_condition = [=](double s, double v, double r) { return std::max(0.0, s - strike); };
//     auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_3d>(terminal_condition);
//     // heat data config:
//     auto const heat_data_ptr = std::make_shared<heat_data_config_3d>(heat_coeffs_data_ptr, heat_init_data_ptr);
//     // spot boundary conditions:
//     auto const &dirichlet_s = [=](double t, double v, double r) { return 0.0; };
//     auto const &neumann_s = [=](double t, double v, double r) { return -1.0; };
//     auto const &boundary_low_s_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_s);
//     auto const &boundary_high_s_ptr = std::make_shared<neumann_boundary_3d>(neumann_s);
//     auto const &x_boundary_pair = std::make_pair(boundary_low_s_ptr, boundary_high_s_ptr);
//     //  upper vol boundary:
//     auto const &dirichlet_v = [=](double t, double s, double r) { return s; };
//     auto const &y_upper_boundary_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_v);
//     // rate boundary conditions:
//     auto const &neumann_low_r = [=](double t, double s, double v) { return 0.0; };
//     auto const &neumann_high_r = [=](double t, double s, double v) { return 0.0; };
//     auto const &boundary_low_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_low_r);
//     auto const &boundary_high_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_high_r);
//     auto const &z_boundary_pair = std::make_pair(boundary_low_r_ptr, boundary_high_r_ptr);
//     // splitting method configuration:
//     auto const &splitting_config_ptr =
//         std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);
//
//     // grid config:
//     auto const alpha_scale = 3.;
//     auto const beta_scale = 50.;
//     auto const &grid_config_hints_ptr =
//         std::make_shared<grid_config_hints_3d>(strike /*, alpha_scale, beta_scale, beta_scale,
//         grid_enum::Nonuniform*/);
//
//     // initialize pde solver
//     hhw_equation pdesolver(heat_data_ptr, discretization_ptr, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair,
//                            splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_0_66_cn_solver_config_ptr);
//     // prepare container for solution:
//     container_3d<by_enum::RowPlane> solution(Sd, Vd, Rd, double{});
//     // get the solution:
//     pdesolver.solve(solution);
//
//     std::stringstream ssn;
//     ssn << "outputs/xmls/" << test_name_num << ".xml";
//     std::string file_name_num{ssn.str()};
//     std::ofstream numerical(file_name_num);
//     xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
//     numerical.close();
//     std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
// }
//
// void test_impl_hhw_dsssolver_dr_cn_srf_xml()
//{
//     std::cout << "============================================================\n";
//     std::cout << "== Implicit Heston-Hull-White Equation (DSS Solver => DR) ==\n";
//     std::cout << "============================================================\n";
//
//     // impl_hhw_dsssolver_dr_cn_srf_xml();
//     impl_hhw_dsssolver_dr_0_66_cn_srf_xml();
//
//     std::cout << "============================================================\n";
// }

////////////////////////////////////////////////////////////////////////////////////
///
///                             EXPLICIT SOLVERS
///
////////////////////////////////////////////////////////////////////////////////////

void expl_heston_host_euler_srf_xml()
{

    const std::string test_name_num{"expl_heston_host_euler_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_expl_bwd_euler_solver_config_ptr;
    using lss_pde_solvers::two_dimensional::explicit_solvers::heston_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Euler Explicit algo\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.8;
    // number of space subdivisions for spot:
    std::size_t const Sd = 50;
    // number of space subdivision for volatility:
    std::size_t const Vd = 30;
    // number of time subdivisions:
    std::size_t const Td = 8000;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              grid_config_hints_ptr, host_expl_bwd_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_expl_heston_host_euler_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "========== Explicit Heston Equation (Euler on Host) =======\n";
    std::cout << "============================================================\n";

    expl_heston_host_euler_srf_xml();

    std::cout << "============================================================\n";
}

void expl_heston_device_euler_srf_xml()
{
    const std::string test_name_num{"expl_heston_device_euler_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::dev_expl_bwd_euler_solver_config_ptr;
    using lss_pde_solvers::two_dimensional::explicit_solvers::heston_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston Call equation: \n\n";
    std::cout << " Using Euler Explicit algo\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,v,t) = 0.5*v*s*s*U_ss(s,v,t) + 0.5*sig*sig*v*U_vv(s,v,t)"
                 " + rho*sig*v*s*U_sv(s,v,t) + r*s*U_s(s,v,t)"
                 " + [k*(theta-v)-lambda*v]*U_v(s,v,t) - r*U(s,v,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(50,v,t) = 0 and  U_s(200,v,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t)+k*theta*U_v(s,0,t)-rU(s,0,t)-U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,v,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.3;
    auto const &sig_kappa = 2.0;
    auto const &sig_theta = 0.2;
    auto const &rho = 0.8;
    // number of space subdivisions for spot:
    std::size_t const Sd = 50;
    // number of space subdivision for volatility:
    std::size_t const Vd = 30;
    // number of time subdivisions:
    std::size_t const Td = 8000;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto a = [=](double t, double s, double v) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v) { return (0.5 * sig_sig * sig_sig * v); };
    auto c = [=](double t, double s, double v) { return (rho * sig_sig * v * s); };
    auto d = [=](double t, double s, double v) { return (rate * s); };
    auto e = [=](double t, double s, double v) { return (sig_kappa * (sig_theta - v)); };
    auto f = [=](double t, double s, double v) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              grid_config_hints_ptr, dev_expl_bwd_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_expl_heston_device_euler_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "========= Explicit Heston Equation (Euler on Device) =======\n";
    std::cout << "============================================================\n";

    expl_heston_device_euler_srf_xml();

    std::cout << "============================================================\n";
}

void expl_sabr_host_euler_srf_xml()
{
    const std::string test_name_num{"expl_sabr_host_euler_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_expl_bwd_euler_solver_config_ptr;
    using lss_pde_solvers::two_dimensional::explicit_solvers::heston_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value SABR Call equation: \n\n";
    std::cout << " Using Double Sweep algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,a,t) = 0.5*a*a*s^(2b)*D^(2*(1-b))*U_ss(s,a,t) "
                 " + 0.5*sig*sig*a*a*U_vv(s,v,t)"
                 " + rho*sig*s^b*D^(1-b)*a*a*U_sv(s,a,t) + r*s*U_s(s,a,t)"
                 " - r*U(s,a,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(0,a,t) = 0 and  U_s(200,a,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t) - rU(s,0,t) - U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,a,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.081;
    auto const &rho = 0.6;
    auto const &beta = 0.7;
    // number of space subdivisions for spot:
    std::size_t const Sd = 50;
    // number of space subdivision for volatility:
    std::size_t const Vd = 40;
    // number of time subdivisions:
    std::size_t const Td = 8000;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto D = [=](double t, double s, double alpha) { return std::exp(-rate * (maturity - t)); };
    auto a = [=](double t, double s, double alpha) {
        return (0.5 * alpha * alpha * std::pow(s, 2.0 * beta) * std::pow(D(t, s, alpha), 2.0 * (1.0 - beta)));
    };
    auto b = [=](double t, double s, double alpha) { return (0.5 * sig_sig * sig_sig * alpha * alpha); };
    auto c = [=](double t, double s, double alpha) {
        return (rho * sig_sig * alpha * alpha * std::pow(s, beta) * std::pow(D(t, s, alpha), (1.0 - beta)));
    };
    auto d = [=](double t, double s, double alpha) { return (rate * s); };
    auto e = [=](double t, double s, double alpha) { return 0.0; };
    auto f = [=](double t, double s, double alpha) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              grid_config_hints_ptr, host_expl_bwd_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_expl_sabr_host_euler_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "=========== Explicit SABR Equation (Euler on Host) =========\n";
    std::cout << "============================================================\n";

    expl_sabr_host_euler_srf_xml();

    std::cout << "============================================================\n";
}

void expl_sabr_device_euler_srf_xml()
{
    const std::string test_name_num{"expl_sabr_device_euler_srf_numerical"};
    using lss_pde_solvers::default_heat_solver_configs::host_expl_bwd_euler_solver_config_ptr;
    using lss_pde_solvers::two_dimensional::explicit_solvers::heston_equation;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value SABR Call equation: \n\n";
    std::cout << " Using Double Sweep algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(s,a,t) = 0.5*a*a*s^(2b)*D^(2*(1-b))*U_ss(s,a,t) "
                 " + 0.5*sig*sig*a*a*U_vv(s,v,t)"
                 " + rho*sig*s^b*D^(1-b)*a*a*U_sv(s,a,t) + r*s*U_s(s,a,t)"
                 " - r*U(s,a,t)\n\n";
    std::cout << " where\n\n";
    std::cout << " 50 < s < 200, 0 < v < 1, and 0 < t < 1,\n";
    std::cout << " U(0,a,t) = 0 and  U_s(200,a,t) - 1 = 0, 0 < t < 1\n";
    std::cout << " r*s*U_s(s,0,t) - rU(s,0,t) - U_t(s,0,t) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(s,1,t) = s, 0 < t < 1\n";
    std::cout << " U(s,a,T) = max(0,s - K), s in <50,200> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &rate = 0.03;
    auto const &sig_sig = 0.081;
    auto const &rho = 0.6;
    auto const &beta = 0.7;
    // number of space subdivisions for spot:
    std::size_t const Sd = 50;
    // number of space subdivision for volatility:
    std::size_t const Vd = 40;
    // number of time subdivisions:
    std::size_t const Td = 8000;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.2);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr =
        std::make_shared<pde_discretization_config_2d>(spacex_range, spacey_range, Sd, Vd, time_range, Td);
    // coeffs:
    auto D = [=](double t, double s, double alpha) { return std::exp(-rate * (maturity - t)); };
    auto a = [=](double t, double s, double alpha) {
        return (0.5 * alpha * alpha * std::pow(s, 2.0 * beta) * std::pow(D(t, s, alpha), 2.0 * (1.0 - beta)));
    };
    auto b = [=](double t, double s, double alpha) { return (0.5 * sig_sig * sig_sig * alpha * alpha); };
    auto c = [=](double t, double s, double alpha) {
        return (rho * sig_sig * alpha * alpha * std::pow(s, beta) * std::pow(D(t, s, alpha), (1.0 - beta)));
    };
    auto d = [=](double t, double s, double alpha) { return (rate * s); };
    auto e = [=](double t, double s, double alpha) { return 0.0; };
    auto f = [=](double t, double s, double alpha) { return (-rate); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_2d>(a, b, c, d, e, f);
    // terminal condition:
    auto terminal_condition = [=](double s, double v) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_2d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_2d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // horizontal spot boundary conditions:
    auto const &dirichlet_low = [=](double t, double v) { return 0.0; };
    auto const &neumann_high = [=](double t, double s) { return -1.0; };
    auto const &boundary_low_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_low);
    auto const &boundary_high_ptr = std::make_shared<neumann_boundary_2d>(neumann_high);
    auto const &horizontal_boundary_pair = std::make_pair(boundary_low_ptr, boundary_high_ptr);
    // vertical upper vol boundary:
    auto const &dirichlet_high = [=](double t, double s) { return s; };
    auto const &vertical_upper_boundary_ptr = std::make_shared<dirichlet_boundary_2d>(dirichlet_high);
    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_2d>(strike, alpha_scale, beta_scale, grid_enum::Nonuniform);

    // initialize pde solver
    heston_equation pdesolver(heat_data_ptr, discretization_ptr, vertical_upper_boundary_ptr, horizontal_boundary_pair,
                              grid_config_hints_ptr, host_expl_bwd_euler_solver_config_ptr);
    // prepare container for solution:
    matrix_2d solution(Sd, Vd, double{});
    // get the solution:
    pdesolver.solve(solution);

    std::stringstream ssn;
    ssn << "outputs/xmls/" << test_name_num << ".xml";
    std::string file_name_num{ssn.str()};
    std::ofstream numerical(file_name_num);
    xml(discretization_ptr, grid_config_hints_ptr, solution, numerical);
    numerical.close();
    std::cout << "Numerical solution has been saved to XML file: " << file_name_num << "\n";
}

void test_expl_sabr_device_euler_srf_xml()
{
    std::cout << "============================================================\n";
    std::cout << "========= Explicit SABR Equation (Euler on Device) =========\n";
    std::cout << "============================================================\n";

    expl_sabr_device_euler_srf_xml();

    std::cout << "============================================================\n";
}
