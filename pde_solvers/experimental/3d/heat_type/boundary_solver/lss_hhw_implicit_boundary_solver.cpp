#include "lss_hhw_implicit_boundary_solver.hpp"

#include "../../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../../discretization/lss_discretization.hpp"
#include "../../../../../discretization/lss_grid.hpp"
#include "../../../../lss_pde_discretization_config.hpp"

#include <iostream>
#include <string>

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_2d;
using lss_boundary::dirichlet_boundary_3d;
using lss_boundary::neumann_boundary_2d;
using lss_boundary::neumann_boundary_3d;
using d_2d = discretization_2d;
using lss_grids::grid_2d;

namespace three_dimensional
{

void implicit_hhw_boundary_rhs::rhs(heat_coefficients_3d_ptr const &cfg, grid_config_3d_ptr const &grid_cfg,
                                    std::size_t const &y_index, double const &y,
                                    boundary_3d_pair const &x_boundary_pair, boundary_3d_pair const &z_boundary_pair,
                                    matrix_3d const &input, double const &time, matrix_2d &solution)
{
    auto const one = 1.0;
    auto const two = 2.0;
    auto const three = 3.0;
    auto const four = 4.0;
    auto const &G = cfg->G_;
    auto const &C = cfg->C_;
    auto const &I = cfg->I_;
    auto const &J = cfg->J_;
    auto const &H = cfg->H_;
    auto const &M_3 = cfg->M_3_;
    auto const &P_3 = cfg->P_3_;
    auto const &theta = cfg->theta_;
    auto const alpha_3 = cfg->alpha_3_;
    auto const gamma_1 = cfg->gamma_1_;
    auto const gamma_2 = cfg->gamma_2_;
    auto const gamma_3 = cfg->gamma_3_;

    auto const &first_x_bnd = x_boundary_pair.first;
    auto const &second_x_bnd = x_boundary_pair.second;
    auto const &first_z_bnd = z_boundary_pair.first;
    auto const &second_z_bnd = z_boundary_pair.second;
    const std::size_t N = solution.rows() - 1;
    const std::size_t L = solution.columns() - 1;

    double x{}, z{}, h_1{}, h_3{};
    h_1 = grid_3d::step_1(grid_cfg);
    h_3 = grid_3d::step_3(grid_cfg);
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(first_x_bnd))
    {
        auto fun = [&](double z) { return ptr->value(time, y, z); };
        container_t lower(double{0.0}, L + 1);
        for (std::size_t k = 0; k < L + 1; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            lower[k] = fun(z);
        }
        solution.row(0) = lower;
    }

    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_x_bnd))
    {
        x = grid_3d::value_1(grid_cfg, N);
        auto fun = [&](double z) { return two * h_1 * ptr->value(time, y, z); };
        container_t upper(double{0.0}, L + 1);
        for (std::size_t k = 1; k < L; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            upper[k] = ((one - theta) * M_3(time, x, y, z) * input(N, y_index, k - 1)) +
                       ((one - theta) * P_3(time, x, y, z) * input(N, y_index, k + 1)) +
                       ((one - three * gamma_3 * H(time, x, y, z) - two * alpha_3 * (one - theta) * C(time, x, y, z)) *
                        input(N, y_index, k)) -
                       (gamma_1 * fun(z) * G(time, x, y, z)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, k)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, k));
        }

        if (auto const &ptr_z = std::dynamic_pointer_cast<neumann_boundary_3d>(first_z_bnd))
        {
            z = grid_3d::value_3(grid_cfg, 0);
            auto const beta_z = two * h_3 * ptr_z->value(time, y, z);

            upper[0] = ((one - theta) * M_3(time, x, y, z) * beta_z) +
                       (two * alpha_3 * (one - theta) * C(time, x, y, z) * input(N, y_index, 1)) +
                       ((one - three * gamma_3 * H(time, x, y, z) - two * alpha_3 * (one - theta) * C(time, x, y, z)) *
                        input(N, y_index, 0)) -
                       (gamma_1 * fun(z) * G(time, x, y, z)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, 0)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, 0));
        }

        if (auto const &ptr_z = std::dynamic_pointer_cast<neumann_boundary_3d>(second_z_bnd))
        {
            z = grid_3d::value_3(grid_cfg, L);
            auto const beta_z = two * h_3 * ptr_z->value(time, y, z);

            upper[L] = (two * alpha_3 * (one - theta) * C(time, x, y, z) * input(N, y_index, L - 1)) -
                       ((one - theta) * P_3(time, x, y, z) * beta_z) +
                       ((one - three * gamma_3 * H(time, x, y, z) - two * alpha_3 * (one - theta) * C(time, x, y, z)) *
                        input(N, y_index, L)) -
                       (gamma_1 * fun(z) * G(time, x, y, z)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, L)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, L));
        }
        solution.row(N) = upper;
    }

    double val{};
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(first_z_bnd))
    {
        z = grid_3d::value_3(grid_cfg, 0);
        auto fun = [&](double x) { return two * h_3 * ptr->value(time, x, y); };
        for (std::size_t i = 1; i < N; ++i)
        {
            x = grid_3d::value_1(grid_cfg, i);
            val = ((one - theta) * fun(x) * M_3(time, x, y, z)) +
                  (two * alpha_3 * (one - theta) * C(time, x, y, z) * input(i, y_index, 1)) +
                  ((one - three * gamma_3 * H(time, x, y, z) - two * alpha_3 * (one - theta) * C(time, x, y, z)) *
                   input(i, y_index, 0)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, 0)) -
                  (gamma_1 * G(time, x, y, z) * input(i - 1, y_index, 0)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, 0)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, 0));
            solution(i, 0, val);
        }
    }

    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_z_bnd))
    {
        z = grid_3d::value_3(grid_cfg, L);
        auto fun = [&](double x) { return two * h_3 * ptr->value(time, x, y); };
        for (std::size_t i = 1; i < N; ++i)
        {
            x = grid_3d::value_1(grid_cfg, i);
            val = (two * alpha_3 * (one - theta) * C(time, x, y, z) * input(i, y_index, L - 1)) -
                  ((one - theta) * fun(x) * P_3(time, x, y, z)) +
                  ((one - three * gamma_3 * H(time, x, y, z) - two * alpha_3 * (one - theta) * C(time, x, y, z)) *
                   input(i, y_index, L)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, L)) -
                  (gamma_1 * G(time, x, y, z) * input(i - 1, y_index, L)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, L)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, L));
            solution(i, L, val);
        }
    }

    for (std::size_t i = 1; i < N; ++i)
    {
        x = grid_3d::value_1(grid_cfg, i);
        for (std::size_t k = 1; k < L; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            val = ((one - theta) * M_3(time, x, y, z) * input(i, y_index, k - 1)) +
                  ((one - theta) * P_3(time, x, y, z) * input(i, y_index, k + 1)) +
                  ((one - three * gamma_3 * H(time, x, y, z) - two * alpha_3 * (one - theta) * C(time, x, y, z)) *
                   input(i, y_index, k)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, k)) -
                  (gamma_1 * G(time, x, y, z) * input(i - 1, y_index, k)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, k)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, k));
            solution(i, k, val);
        }
    }
}

void implicit_hhw_boundary_rhs::rhs_source(heat_coefficients_3d_ptr const &cfg, grid_config_3d_ptr const &grid_cfg,
                                           std::size_t const &y_index, double const &y,
                                           boundary_3d_pair const &x_boundary_pair,
                                           boundary_3d_pair const &z_boundary_pair, matrix_3d const &input,
                                           matrix_3d const &inhom_input, double const &time, matrix_2d &solution)
{
    auto const one = 1.0;
    auto const two = 2.0;
    auto const three = 3.0;
    auto const four = 4.0;
    auto const &G = cfg->G_;
    auto const &C = cfg->C_;
    auto const &I = cfg->I_;
    auto const &J = cfg->J_;
    auto const &H = cfg->H_;
    auto const &M_3 = cfg->M_3_;
    auto const &P_3 = cfg->P_3_;
    auto const alpha_1 = cfg->alpha_1_;
    auto const alpha_3 = cfg->alpha_3_;
    auto const gamma_1 = cfg->gamma_1_;
    auto const gamma_2 = cfg->gamma_2_;
    auto const gamma_3 = cfg->gamma_3_;
    auto const rho = cfg->rho_;

    auto const &first_x_bnd = x_boundary_pair.first;
    auto const &second_x_bnd = x_boundary_pair.second;
    auto const &first_z_bnd = z_boundary_pair.first;
    auto const &second_z_bnd = z_boundary_pair.second;
    const std::size_t N = solution.rows() - 1;
    const std::size_t L = solution.columns() - 1;

    double x{}, z{}, h_1{}, h_3{};
    h_1 = grid_3d::step_1(grid_cfg);
    h_3 = grid_3d::step_3(grid_cfg);
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(first_x_bnd))
    {
        auto fun = [&](double z) { return ptr->value(time, y, z); };
        container_t lower(double{0.0}, L);
        for (std::size_t k = 0; k < L + 1; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            lower[k] = fun(z);
        }
        solution.row(0) = lower;
    }
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_x_bnd))
    {
        x = grid_3d::value_1(grid_cfg, N);
        auto fun = [&](double z) { return two * h_1 * ptr->value(time, y, z); };
        container_t upper(double{0.0}, L);
        for (std::size_t k = 0; k < L; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            upper[k] = (-gamma_1 * fun(z) * G(time, x, y, z)) + (M_3(time, x, y, z) * input(N, y_index, k - 1)) +
                       ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_3 * H(time, x, y, z) +
                         rho * J(time, x, y, z)) *
                        input(N, y_index, k)) +
                       (P_3(time, x, y, z) * input(N, y_index, k + 1)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, k)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, k)) + rho * inhom_input(N, y_index, k);
        }

        if (auto const &ptr_z = std::dynamic_pointer_cast<neumann_boundary_3d>(first_z_bnd))
        {
            z = grid_3d::value_3(grid_cfg, 0);
            auto const beta_z = two * h_3 * ptr_z->value(time, y, z);

            upper[0] = (-gamma_1 * fun(z) * G(time, x, y, z)) + (M_3(time, x, y, z) * beta_z) +
                       ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_3 * H(time, x, y, z) +
                         rho * J(time, x, y, z)) *
                        input(N, y_index, 0)) +
                       (two * alpha_3 * C(time, x, y, z) * input(N, y_index, 1)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, 0)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, 0)) + rho * inhom_input(N, y_index, 0);
        }

        if (auto const &ptr_z = std::dynamic_pointer_cast<neumann_boundary_3d>(second_z_bnd))
        {
            z = grid_3d::value_3(grid_cfg, L);
            auto const beta_z = two * h_3 * ptr_z->value(time, y, z);

            upper[L] = (-gamma_1 * fun(z) * G(time, x, y, z)) - (P_3(time, x, y, z) * beta_z) +
                       ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_3 * H(time, x, y, z) +
                         rho * J(time, x, y, z)) *
                        input(N, y_index, L)) +
                       (two * alpha_3 * C(time, x, y, z) * input(N, y_index, L - 1)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, L)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, L)) + rho * inhom_input(N, y_index, L);
        }
        solution.row(N) = upper;
    }

    double val{};
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(first_z_bnd))
    {
        z = grid_3d::value_3(grid_cfg, 0);
        auto fun = [&](double x) { return two * h_3 * ptr->value(time, x, y); };
        for (std::size_t i = 1; i < N; ++i)
        {
            x = grid_3d::value_1(grid_cfg, i);
            val = (-gamma_1 * G(time, x, y, z) * input(i - 1, y_index, 0)) + (M_3(time, x, y, z) * fun(x)) +
                  ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_2 * H(time, x, y, z) +
                    rho * J(time, x, y, z)) *
                   input(i, y_index, 0)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, 0)) +
                  (two * alpha_3 * C(time, x, y, z) * input(i, y_index, 1)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, 0)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, 0)) + rho * inhom_input(i, y_index, 0);
            solution(i, 0, val);
        }
    }

    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_z_bnd))
    {
        z = grid_3d::value_3(grid_cfg, L);
        auto fun = [&](double x) { return two * h_3 * ptr->value(time, x, y); };
        for (std::size_t i = 1; i < N; ++i)
        {
            x = grid_3d::value_1(grid_cfg, i);
            val = (-gamma_1 * G(time, x, y, z) * input(i - 1, y_index, L)) - (P_3(time, x, y, z) * fun(x)) +
                  ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_2 * H(time, x, y, z) +
                    rho * J(time, x, y, z)) *
                   input(i, y_index, L)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, L)) +
                  (two * alpha_3 * C(time, x, y, z) * input(i, y_index, L - 1)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, L)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, L)) + rho * inhom_input(i, y_index, L);
            solution(i, L, val);
        }
    }

    for (std::size_t i = 1; i < N; ++i)
    {
        x = grid_3d::value_1(grid_cfg, i);
        for (std::size_t k = 1; k < L; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            val = (-gamma_1 * G(time, x, y, z) * input(i - 1, y_index, k)) +
                  (M_3(time, x, y, z) * input(i, y_index, k - 1)) +
                  ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_2 * H(time, x, y, z) +
                    rho * J(time, x, y, z)) *
                   input(i, y_index, k)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, k)) +
                  (P_3(time, x, y, z) * input(i, y_index, k + 1)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, k)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, k)) + rho * inhom_input(i, y_index, k);
            solution(i, k, val);
        }
    }
}

void hhw_implicit_boundary_solver::split(double const &x, double const &y, double const &time, container_t &low,
                                         container_t &diag, container_t &high)
{
    double z{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        z = grid_3d::value_3(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->M_3_(time, x, y, z));
        diag[t] = (cone_ + ctwo_ * coefficients_->theta_ * coefficients_->alpha_3_ * coefficients_->C_(time, x, y, z));
        high[t] = (-coefficients_->theta_ * coefficients_->P_3_(time, x, y, z));
    }
}

hhw_implicit_boundary_solver::hhw_implicit_boundary_solver(
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru, heat_coefficients_3d_ptr const &coefficients,
    grid_config_3d_ptr const &grid_config)
    : solveru_ptr_{solveru}, coefficients_{coefficients}, grid_cfg_{grid_config}
{
}

hhw_implicit_boundary_solver::~hhw_implicit_boundary_solver()
{
}

void hhw_implicit_boundary_solver::solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
                                         boundary_3d_ptr const &y_upper_boundary_ptr,
                                         boundary_3d_pair const &z_boundary_pair, double const &time,
                                         matrix_3d &solution)
{
    // 2D container for intermediate solution:
    matrix_2d solution_p(coefficients_->space_size_x_, coefficients_->space_size_z_, double{});
    /// get the right-hand side of the scheme:
    auto const y = grid_3d::value_2(grid_cfg_, 0);
    implicit_hhw_boundary_rhs::rhs(coefficients_, grid_cfg_, 0, y, x_boundary_pair, z_boundary_pair, prev_solution,
                                   time, solution_p);
    // 1D container for intermediate solution:
    container_t solution_v(double{}, coefficients_->space_size_z_);
    // containers for second split solver:
    low_.resize(coefficients_->space_size_z_);
    diag_.resize(coefficients_->space_size_z_);
    high_.resize(coefficients_->space_size_z_);
    double x{};
    for (std::size_t i = 1; i < coefficients_->space_size_x_ - 1; ++i)
    {
        x = grid_3d::value_1(grid_cfg_, i);
        split(x, y, time, low_, diag_, high_);
        solveru_ptr_->set_diagonals(low_, diag_, high_);
        solveru_ptr_->set_rhs(solution_p.row(i));
        solveru_ptr_->solve(z_boundary_pair, solution_v, time, x, y);
        solution_p.row(i) = solution_v;
    }
    // boundary plane at Y = 0:
    solution.column_plane(0) = solution_p.data();

    auto const &upper_bnd_ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(y_upper_boundary_ptr);
    auto const &upper_bnd = [=](double t, double s, double r) { return upper_bnd_ptr->value(t, s, r); };
    d_2d::of_function(grid_cfg_->grid_13(), time, upper_bnd, solution_p);
    // boundary plane at Y = M-1:
    solution.column_plane(coefficients_->space_size_y_ - 1) = solution_p.data();
}

void hhw_implicit_boundary_solver::solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
                                         boundary_3d_pair const &z_boundary_pair, double const &time,
                                         matrix_3d &solution)
{
    // 3D container for intermediate solution:
    matrix_2d solution_yz(coefficients_->space_size_y_, coefficients_->space_size_z_, double{});
    matrix_2d solution_xy(coefficients_->space_size_x_, coefficients_->space_size_y_, double{});
    // some constants:
    auto const &start_y = coefficients_->rangey_->lower();
    // prepare grid_xy:
    auto const &grid_12 = grid_cfg_->grid_12();
    // prepare grid_yz:
    auto const &grid_23 = grid_cfg_->grid_23();

    /// populating lower X:
    auto const &lower_x_ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(x_boundary_pair.first);
    auto const &lower_x_bnd = [=](double t, double v, double r) { return lower_x_ptr->value(t, v, r); };
    d_2d::of_function(grid_23, time, lower_x_bnd, solution_yz);
    solution.row_plane(0) = solution_yz.data();

    // populating upper X:
    auto const lri = solution.rows() - 1;
    auto const lci = solution.columns() - 1;
    auto const lli = solution.layers() - 1;
    auto const two = 2.0;
    auto const three = 3.0;
    auto const four = 4.0;
    auto const &upper_x_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(x_boundary_pair.second);
    auto const &upper_x_bnd = [=](double t, double v, double r) {
        const std::size_t j = grid_3d::index_of_2(grid_cfg_, v);
        const std::size_t k = grid_3d::index_of_3(grid_cfg_, r);
        auto const bnd_val = upper_x_ptr->value(t, v, r);
        auto const h_1 = grid_3d::step_1(grid_cfg_);
        return (((four * solution(lri - 1, j, k)) - solution(lri - 2, j, k) - (two * h_1 * bnd_val)) / three);
    };
    d_2d::of_function(grid_23, time, upper_x_bnd, solution_yz);
    solution.row_plane(lri) = solution_yz.data();

    // populating lower Z:
    auto const &lower_z_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(z_boundary_pair.first);
    auto const &lower_z_bnd = [=](double t, double s, double v) {
        const std::size_t i = grid_3d::index_of_1(grid_cfg_, s);
        const std::size_t j = grid_3d::index_of_2(grid_cfg_, v);
        auto const bnd_val = lower_z_ptr->value(t, s, v);
        auto const h_3 = grid_3d::step_3(grid_cfg_);
        return (((four * solution(i, j, 1)) - solution(i, j, 2) + (two * h_3 * bnd_val)) / three);
    };
    d_2d::of_function(grid_12, time, lower_z_bnd, solution_xy);
    solution.layer_plane(0) = solution_xy.data();

    // populating upper Z:
    auto const &upper_z_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(z_boundary_pair.second);
    auto const &upper_z_bnd = [=](double t, double s, double v) {
        const std::size_t i = grid_3d::index_of_1(grid_cfg_, s);
        const std::size_t j = grid_3d::index_of_2(grid_cfg_, v);
        auto const bnd_val = upper_z_ptr->value(t, s, v);
        auto const h_3 = grid_3d::step_3(grid_cfg_);
        return (((four * solution(i, j, lli - 1)) - solution(i, j, lli - 2) - (two * h_3 * bnd_val)) / three);
    };
    d_2d::of_function(grid_12, time, upper_z_bnd, solution_xy);
    solution.layer_plane(lli) = solution_xy.data();
}

} // namespace three_dimensional

} // namespace lss_pde_solvers
