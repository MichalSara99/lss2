#if !defined(_LSS_HHW_EQUATION_T_HPP_)
#define _LSS_HHW_EQUATION_T_HPP_

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_print.hpp"
#include "../../../common/lss_range.hpp"
#include "../../../containers/lss_matrix_3d.hpp"
#include "pde_solvers/experimental/3d/heat_type/lss_hhw_equation.hpp"
#include <map>

using lss_boundary::dirichlet_boundary_3d;
using lss_boundary::neumann_boundary_3d;
using lss_containers::matrix_3d;
using lss_enumerations::factorization_enum;
using lss_enumerations::memory_space_enum;
using lss_enumerations::splitting_method_enum;
using lss_enumerations::traverse_direction_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_grids::grid_config_hints_3d;
using lss_pde_solvers::heat_coefficient_data_config_3d;
using lss_pde_solvers::heat_data_config_3d;
using lss_pde_solvers::heat_implicit_solver_config;
using lss_pde_solvers::heat_initial_data_config_3d;
using lss_pde_solvers::pde_discretization_config_3d;
using lss_pde_solvers::splitting_method_config;
using lss_pde_solvers::three_dimensional::implicit_solvers::hhw_equation;
using lss_print::print;
using lss_utility::range;

// ///////////////////////////////////////////////////////////////////////////
//							HESTON-HULL-WHITE PROBLEMS
// ///////////////////////////////////////////////////////////////////////////

// ===========================================================================
// ========================== IMPLICIT SOLVERS ===============================
// ===========================================================================

void impl_hhw_equation_cuda_qr_solver_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::dev_bwd_cusolver_qr_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston-Hull-White Call equation: \n\n";
    std::cout << "Using DSS algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(t,s,v,r) = 0.5*s*s*v*U_ss(t,s,v,r) +"
                 "0.5*vol_1*vol_1*v*U_vv(t,s,v,r) +"
                 "0.5*vol_2*vol_2*U_rr(t,s,v,r) + "
                 "rho_12*vol_1*s*v*U_sv(t,s,v,r) + "
                 "rho_13*vol_2*s*sqrt(v)*U_sr(t,s,v,r) + "
                 "rho_23*vol_1*vol_2*sqrt(v)*U_vr(t,s,v,r) + "
                 "r*s*U_s(t,s,v,r) + [k*(theta-v)-lambda*v]*U_v(t,s,v,r) +"
                 "a*(b-r)*U_r(t,s,v,r) - r*U(t,s,v,r)\n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < s < 20, 0 < v < 1,-1 < r < 1, and 0 < t < 1,\n";
    std::cout << " U(t,0,v,r) = 0 and  U_s(t,20,v,r) - 1 = 0, 0 < t < 1\n";
    std::cout << " 0.5*vol_2*vol_2*U_rr(t,s,0,r) + r*s*U_s(t,s,0,r) + "
                 "k*theta*U_v(t,s,0,r) + a*(b-r)*U_r(t,s,0,r) - "
                 "rU(t,s,0,r) - U_t(t,s,0,r) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(t,s,1,r) = s, 0 < t < 1\n";
    std::cout << " U_r(t,s,v,-1) = 0 and  U_r(t,s,v,1) = 0, 0 < t < 1\n";
    std::cout << " U(T,s,v,r) = max(0,s - K), s in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &v_sig = 0.8;
    auto const &v_kappa = 3.0;
    auto const &v_theta = 0.2;
    auto const &rho_12 = 0.6;
    auto const &rho_13 = 0.2;
    auto const &rho_23 = 0.4;
    auto const &c_1 = 0.3;
    auto const &c_2 = 0.01;
    auto const &c_3 = 0.02;
    auto const &r_a = 0.2;
    auto const &r_sig = 0.05;
    // number of space subdivisions for spot:
    std::size_t const Sd = 50;
    // number of space subdivision for volatility:
    std::size_t const Vd = 25;
    // number of space subdivision for rate:
    std::size_t const Rd = 25;
    // number of time subdivisions:
    std::size_t const Td = 500;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
    // space Rate range:
    auto const &spacez_range = std::make_shared<range>(0.12, 0.21);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_3d>(
        spacex_range, spacey_range, spacez_range, Sd, Vd, Rd, time_range, Td);
    // coeffs:
    auto B = [=](double t) { return (c_1 - c_2 * std::exp(-c_3 * (maturity - t))); };
    auto a = [=](double t, double s, double v, double r) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v, double r) { return (0.5 * v_sig * v_sig * v); };
    auto c = [=](double t, double s, double v, double r) { return (0.5 * r_sig * r_sig); };
    auto d = [=](double t, double s, double v, double r) { return (rho_12 * v_sig * s * v); };
    auto e = [=](double t, double s, double v, double r) { return (rho_13 * r_sig * s * sqrt(v)); };
    auto f = [=](double t, double s, double v, double r) { return (rho_23 * v_sig * r_sig * sqrt(v)); };
    auto g = [=](double t, double s, double v, double r) { return (r * s); };
    auto h = [=](double t, double s, double v, double r) { return (v_kappa * (v_theta - v)); };
    auto i = [=](double t, double s, double v, double r) { return (r_a * (B(t) - r)); };
    auto j = [=](double t, double s, double v, double r) { return (-r); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_3d>(a, b, c, d, e, f, g, h, i, j);
    // terminal condition:
    auto terminal_condition = [=](double s, double v, double r) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_3d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_3d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // spot boundary conditions:
    auto const &dirichlet_s = [=](double t, double v, double r) { return 0.0; };
    auto const &neumann_s = [=](double t, double v, double r) { return -1.0; };
    auto const &boundary_low_s_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_s);
    auto const &boundary_high_s_ptr = std::make_shared<neumann_boundary_3d>(neumann_s);
    auto const &x_boundary_pair = std::make_pair(boundary_low_s_ptr, boundary_high_s_ptr);
    //  upper vol boundary:
    auto const &dirichlet_v = [=](double t, double s, double r) { return s; };
    auto const &y_upper_boundary_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_v);
    // rate boundary conditions:
    auto const &neumann_low_r = [=](double t, double s, double v) { return 0.0; };
    auto const &neumann_high_r = [=](double t, double s, double v) { return 0.0; };
    auto const &boundary_low_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_low_r);
    auto const &boundary_high_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_high_r);
    auto const &z_boundary_pair = std::make_pair(boundary_low_r_ptr, boundary_high_r_ptr);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);

    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_3d>(strike);

    // initialize pde solver
    hhw_equation pdesolver(heat_data_ptr, discretization_ptr, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair,
                           splitting_config_ptr, grid_config_hints_ptr, dev_bwd_cusolver_qr_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_3d solution(Sd, Vd, Rd, double{});
    // get the solution:
    pdesolver.solve(solution);

    print(discretization_ptr, grid_config_hints_ptr, solution);
}

void test_impl_hhw_equation_cuda_qr_solver()
{
    std::cout << "============================================================\n";
    std::cout << "=============== Implicit HHW Equation (Dir BC) =============\n";
    std::cout << "============================================================\n";

    impl_hhw_equation_cuda_qr_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

void impl_hhw_equation_thomas_lu_solver_crank_nicolson()
{
    using lss_pde_solvers::default_heat_solver_configs::host_bwd_tlusolver_cn_solver_config_ptr;

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston-Hull-White Call equation: \n\n";
    std::cout << "Using Thomas LU algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(t,s,v,r) = 0.5*s*s*v*U_ss(t,s,v,r) +"
                 "0.5*vol_1*vol_1*v*U_vv(t,s,v,r) +"
                 "0.5*vol_2*vol_2*U_rr(t,s,v,r) + "
                 "rho_12*vol_1*s*v*U_sv(t,s,v,r) + "
                 "rho_13*vol_2*s*sqrt(v)*U_sr(t,s,v,r) + "
                 "rho_23*vol_1*vol_2*sqrt(v)*U_vr(t,s,v,r) + "
                 "r*s*U_s(t,s,v,r) + [k*(theta-v)-lambda*v]*U_v(t,s,v,r) +"
                 "a*(b-r)*U_r(t,s,v,r) - r*U(t,s,v,r)\n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < s < 20, 0 < v < 1,-1 < r < 1, and 0 < t < 1,\n";
    std::cout << " U(t,0,v,r) = 0 and  U_s(t,20,v,r) - 1 = 0, 0 < t < 1\n";
    std::cout << " 0.5*vol_2*vol_2*U_rr(t,s,0,r) + r*s*U_s(t,s,0,r) + "
                 "k*theta*U_v(t,s,0,r) + a*(b-r)*U_r(t,s,0,r) - "
                 "rU(t,s,0,r) - U_t(t,s,0,r) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(t,s,1,r) = s, 0 < t < 1\n";
    std::cout << " U_r(t,s,v,-1) = 0 and  U_r(t,s,v,1) = 0, 0 < t < 1\n";
    std::cout << " U(T,s,v,r) = max(0,s - K), s in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 1.0;
    auto const &v_sig = 0.8;
    auto const &v_kappa = 3.0;
    auto const &v_theta = 0.2;
    auto const &rho_12 = 0.6;
    auto const &rho_13 = 0.2;
    auto const &rho_23 = 0.4;
    auto const &c_1 = 0.3;
    auto const &c_2 = 0.01;
    auto const &c_3 = 0.02;
    auto const &r_a = 0.2;
    auto const &r_sig = 0.05;
    // number of space subdivisions for spot:
    std::size_t const Sd = 30;
    // number of space subdivision for volatility:
    std::size_t const Vd = 18;
    // number of space subdivision for rate:
    std::size_t const Rd = 18;
    // number of time subdivisions:
    std::size_t const Td = 500;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.0, 1.0);
    // space Rate range:
    auto const &spacez_range = std::make_shared<range>(0.12, 0.21);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_3d>(
        spacex_range, spacey_range, spacez_range, Sd, Vd, Rd, time_range, Td);
    // coeffs:
    auto B = [=](double t) { return (c_1 - c_2 * std::exp(-c_3 * (maturity - t))); };
    auto a = [=](double t, double s, double v, double r) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v, double r) { return (0.5 * v_sig * v_sig * v); };
    auto c = [=](double t, double s, double v, double r) { return (0.5 * r_sig * r_sig); };
    auto d = [=](double t, double s, double v, double r) { return (rho_12 * v_sig * s * v); };
    auto e = [=](double t, double s, double v, double r) { return (rho_13 * r_sig * s * sqrt(v)); };
    auto f = [=](double t, double s, double v, double r) { return (rho_23 * v_sig * r_sig * sqrt(v)); };
    auto g = [=](double t, double s, double v, double r) { return (r * s); };
    auto h = [=](double t, double s, double v, double r) { return (v_kappa * (v_theta - v)); };
    auto i = [=](double t, double s, double v, double r) { return (r_a * (B(t) - r)); };
    auto j = [=](double t, double s, double v, double r) { return (-r); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_3d>(a, b, c, d, e, f, g, h, i, j);
    // terminal condition:
    auto terminal_condition = [=](double s, double v, double r) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_3d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_3d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // spot boundary conditions:
    auto const &dirichlet_s = [=](double t, double v, double r) { return 0.0; };
    auto const &neumann_s = [=](double t, double v, double r) { return -1.0; };
    auto const &boundary_low_s_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_s);
    auto const &boundary_high_s_ptr = std::make_shared<neumann_boundary_3d>(neumann_s);
    auto const &x_boundary_pair = std::make_pair(boundary_low_s_ptr, boundary_high_s_ptr);
    //  upper vol boundary:
    auto const &dirichlet_v = [=](double t, double s, double r) { return s; };
    auto const &y_upper_boundary_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_v);
    // rate boundary conditions:
    auto const &neumann_low_r = [=](double t, double s, double v) { return 0.0; };
    auto const &neumann_high_r = [=](double t, double s, double v) { return 0.0; };
    auto const &boundary_low_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_low_r);
    auto const &boundary_high_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_high_r);
    auto const &z_boundary_pair = std::make_pair(boundary_low_r_ptr, boundary_high_r_ptr);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);

    // grid config:
    auto const &grid_config_hints_ptr = std::make_shared<grid_config_hints_3d>(strike);

    // initialize pde solver
    hhw_equation pdesolver(heat_data_ptr, discretization_ptr, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair,
                           splitting_config_ptr, grid_config_hints_ptr, host_bwd_tlusolver_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_3d solution(Sd, Vd, Rd, double{});
    // get the solution:
    pdesolver.solve(solution);

    print(discretization_ptr, grid_config_hints_ptr, solution);
}

void impl_hhw_equation_dsssolver_solver_crank_nicolson()
{

    using lss_pde_solvers::default_heat_solver_configs::build_implicit_config;

    auto const host_bwd_dssolver_0_66_cn_solver_config_ptr =
        build_implicit_config(memory_space_enum::Host, traverse_direction_enum::Backward,
                              tridiagonal_method_enum::ThomasLUSolver, factorization_enum::None, (2.0 / 3.0));

    std::cout << "============================================================\n";
    std::cout << "Solving Boundary-value Heston-Hull-White Call equation: \n\n";
    std::cout << "Using DSS algo with implicit Crank-Nicolson method\n\n";
    std::cout << " Value type: " << typeid(double).name() << "\n\n";
    std::cout << " U_t(t,s,v,r) = 0.5*s*s*v*U_ss(t,s,v,r) +"
                 "0.5*vol_1*vol_1*v*U_vv(t,s,v,r) +"
                 "0.5*vol_2*vol_2*U_rr(t,s,v,r) + "
                 "rho_12*vol_1*s*v*U_sv(t,s,v,r) + "
                 "rho_13*vol_2*s*sqrt(v)*U_sr(t,s,v,r) + "
                 "rho_23*vol_1*vol_2*sqrt(v)*U_vr(t,s,v,r) + "
                 "r*s*U_s(t,s,v,r) + [k*(theta-v)-lambda*v]*U_v(t,s,v,r) +"
                 "a*(b-r)*U_r(t,s,v,r) - r*U(t,s,v,r)\n\n";
    std::cout << " where\n\n";
    std::cout << " 0 < s < 20, 0 < v < 1,-1 < r < 1, and 0 < t < 1,\n";
    std::cout << " U(t,0,v,r) = 0 and  U_s(t,20,v,r) - 1 = 0, 0 < t < 1\n";
    std::cout << " 0.5*vol_2*vol_2*U_rr(t,s,0,r) + r*s*U_s(t,s,0,r) + "
                 "k*theta*U_v(t,s,0,r) + a*(b-r)*U_r(t,s,0,r) - "
                 "rU(t,s,0,r) - U_t(t,s,0,r) = 0,"
                 "0 < t < 1\n";
    std::cout << " U(t,s,1,r) = s, 0 < t < 1\n";
    std::cout << " U_r(t,s,v,-1) = 0 and  U_r(t,s,v,1) = 0, 0 < t < 1\n";
    std::cout << " U(T,s,v,r) = max(0,s - K), s in <0,20> \n\n";
    std::cout << "============================================================\n";

    // set up call option parameters:
    auto const &strike = 100.0;
    auto const &maturity = 0.25;
    auto const &v_sig = 0.5;
    auto const &v_kappa = 2.5;
    auto const &v_theta = 0.06;
    auto const &rho_12 = -0.1;
    auto const &rho_13 = -0.3;
    auto const &rho_23 = 0.2;
    auto const &c_1 = 0.005;
    auto const &c_2 = 0.001;
    auto const &c_3 = 2.3;
    auto const &r_a = 0.15;
    auto const &r_sig = 0.1;

    // number of space subdivisions for spot:
    std::size_t const Sd = 75;
    // number of space subdivision for volatility:
    std::size_t const Vd = 35;
    // number of space subdivision for rate:
    std::size_t const Rd = 35;
    // number of time subdivisions:
    std::size_t const Td = 600;
    // space Spot range:
    auto const &spacex_range = std::make_shared<range>(50.0, 200.0);
    // space Vol range:
    auto const &spacey_range = std::make_shared<range>(0.05, 1.0);
    // space Rate range:
    auto const &spacez_range = std::make_shared<range>(0.001, 0.01);
    // time range
    auto const &time_range = std::make_shared<range>(0.0, maturity);
    // discretization config:
    auto const discretization_ptr = std::make_shared<pde_discretization_config_3d>(
        spacex_range, spacey_range, spacez_range, Sd, Vd, Rd, time_range, Td);
    // coeffs:
    auto B = [=](double t) { return (c_1 - c_2 * std::exp(-c_3 * (maturity - t))); };
    auto a = [=](double t, double s, double v, double r) { return (0.5 * v * s * s); };
    auto b = [=](double t, double s, double v, double r) { return (0.5 * v_sig * v_sig * v); };
    auto c = [=](double t, double s, double v, double r) { return (0.5 * r_sig * r_sig); };
    auto d = [=](double t, double s, double v, double r) { return (rho_12 * v_sig * s * v); };
    auto e = [=](double t, double s, double v, double r) { return (rho_13 * r_sig * s * sqrt(v)); };
    auto f = [=](double t, double s, double v, double r) { return (rho_23 * v_sig * r_sig * sqrt(v)); };
    auto g = [=](double t, double s, double v, double r) { return (r * s); };
    auto h = [=](double t, double s, double v, double r) { return (v_kappa * (v_theta - v)); };
    auto i = [=](double t, double s, double v, double r) { return (r_a * (B(t) - r)); };
    auto j = [=](double t, double s, double v, double r) { return (-r); };
    auto const heat_coeffs_data_ptr = std::make_shared<heat_coefficient_data_config_3d>(a, b, c, d, e, f, g, h, i, j);
    // terminal condition:
    auto terminal_condition = [=](double s, double v, double r) { return std::max(0.0, s - strike); };
    auto const heat_init_data_ptr = std::make_shared<heat_initial_data_config_3d>(terminal_condition);
    // heat data config:
    auto const heat_data_ptr = std::make_shared<heat_data_config_3d>(heat_coeffs_data_ptr, heat_init_data_ptr);
    // spot boundary conditions:
    auto const &dirichlet_s = [=](double t, double v, double r) { return 0.0; };
    auto const &neumann_s = [=](double t, double v, double r) { return -1.0; };
    auto const &boundary_low_s_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_s);
    auto const &boundary_high_s_ptr = std::make_shared<neumann_boundary_3d>(neumann_s);
    auto const &x_boundary_pair = std::make_pair(boundary_low_s_ptr, boundary_high_s_ptr);
    //  upper vol boundary:
    auto const &dirichlet_v = [=](double t, double s, double r) { return s; };
    auto const &y_upper_boundary_ptr = std::make_shared<dirichlet_boundary_3d>(dirichlet_v);
    // rate boundary conditions:
    auto const &neumann_low_r = [=](double t, double s, double v) { return 0.0; };
    auto const &neumann_high_r = [=](double t, double s, double v) { return 0.0; };
    auto const &boundary_low_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_low_r);
    auto const &boundary_high_r_ptr = std::make_shared<neumann_boundary_3d>(neumann_high_r);
    auto const &z_boundary_pair = std::make_pair(boundary_low_r_ptr, boundary_high_r_ptr);
    // splitting method configuration:
    auto const &splitting_config_ptr =
        std::make_shared<splitting_method_config>(splitting_method_enum::DouglasRachford);

    // grid config:
    auto const alpha_scale = 3.;
    auto const beta_scale = 50.;
    auto const &grid_config_hints_ptr =
        std::make_shared<grid_config_hints_3d>(strike /*, alpha_scale, beta_scale, beta_scale, grid_enum::Nonuniform*/);

    // initialize pde solver
    hhw_equation pdesolver(heat_data_ptr, discretization_ptr, x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair,
                           splitting_config_ptr, grid_config_hints_ptr, host_bwd_dssolver_0_66_cn_solver_config_ptr);
    // prepare container for solution:
    matrix_3d solution(Sd, Vd, Rd, double{});
    // get the solution:
    pdesolver.solve(solution);

    print(discretization_ptr, grid_config_hints_ptr, solution);
}

void test_impl_hhw_equation_tlu_dss_solver()
{
    std::cout << "============================================================\n";
    std::cout << "=============== Implicit HHW Equation (Dir BC) =============\n";
    std::cout << "============================================================\n";

    // impl_hhw_equation_thomas_lu_solver_crank_nicolson();
    impl_hhw_equation_dsssolver_solver_crank_nicolson();

    std::cout << "============================================================\n";
}

#endif //_LSS_HHW_EQUATION_T_HPP_
