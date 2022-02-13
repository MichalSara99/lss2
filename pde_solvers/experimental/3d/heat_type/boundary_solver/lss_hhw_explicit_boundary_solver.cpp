#include "lss_hhw_explicit_boundary_solver.hpp"

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
using d_2d = discretization_2d<std::vector, std::allocator<double>>;
using lss_grids::grid_2d;

namespace three_dimensional
{

void explicit_hhw_boundary_scheme::rhs(hhw_implicit_coefficients_ptr const &cfg, grid_config_3d_ptr const &grid_cfg,
                                       std::size_t const &y_index, double const &y,
                                       boundary_3d_pair const &x_boundary_pair, boundary_3d_pair const &z_boundary_pair,
                                       container_3d<by_enum::RowPlane> const &input, double const &time,
                                       container_2d<by_enum::Row> &solution)
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
        container_t lower(L + 1, double{0.0});
        for (std::size_t k = 0; k < L + 1; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            lower[k] = fun(z);
        }
        solution(0, lower);
    }

    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_x_bnd))
    {
        x = grid_3d::value_1(grid_cfg, N);
        auto fun = [&](double z) { return two * h_1 * ptr->value(time, y, z); };
        container_t upper(L + 1, double{0.0});
        for (std::size_t k = 1; k < L; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            upper[k] = (-gamma_1 * fun(z) * G(time, x, y, z)) - (M_3(time, x, y, z) * input(N, y_index, k - 1)) +
                       ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_3 * H(time, x, y, z) +
                         rho * J(time, x, y, z)) *
                        input(N, y_index, k)) +
                       (P_3(time, x, y, z) * input(N, y_index, k + 1)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, k)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, k));
        }

        if (auto const &ptr_z = std::dynamic_pointer_cast<neumann_boundary_3d>(first_z_bnd))
        {
            z = grid_3d::value_3(grid_cfg, 0);
            auto const beta_z = two * h_3 * ptr_z->value(time, y, z);

            upper[0] = (-gamma_1 * fun(z) * G(time, x, y, z)) - (M_3(time, x, y, z) * beta_z) +
                       ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_3 * H(time, x, y, z) +
                         rho * J(time, x, y, z)) *
                        input(N, y_index, 0)) +
                       (two * alpha_3 * C(time, x, y, z) * input(N, y_index, 1)) +
                       (four * gamma_2 * H(time, x, y, z) * input(N, y_index + 1, 0)) -
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, 0));
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
                       (gamma_2 * H(time, x, y, z) * input(N, y_index + 2, L));
        }
        solution(N, upper);
    }

    double val{};
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(first_z_bnd))
    {
        z = grid_3d::value_3(grid_cfg, 0);
        auto fun = [&](double x) { return two * h_3 * ptr->value(time, x, y); };
        for (std::size_t i = 1; i < N; ++i)
        {
            x = grid_3d::value_1(grid_cfg, i);
            val = (-gamma_1 * G(time, x, y, z) * input(i - 1, y_index, 0)) - (M_3(time, x, y, z) * fun(x)) +
                  ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_2 * H(time, x, y, z) +
                    rho * J(time, x, y, z)) *
                   input(i, y_index, 0)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, 0)) +
                  (two * alpha_3 * C(time, x, y, z) * input(i, y_index, 1)) +
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
            val = (-gamma_1 * G(time, x, y, z) * input(i - 1, y_index, L)) - (P_3(time, x, y, z) * fun(x)) +
                  ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_2 * H(time, x, y, z) +
                    rho * J(time, x, y, z)) *
                   input(i, y_index, L)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, L)) +
                  (two * alpha_3 * C(time, x, y, z) * input(i, y_index, L - 1)) +
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
            val = (-gamma_1 * G(time, x, y, z) * input(i - 1, y_index, k)) -
                  (M_3(time, x, y, z) * input(i, y_index, k - 1)) +
                  ((one - two * alpha_3 * C(time, x, y, z) - three * gamma_2 * H(time, x, y, z) +
                    rho * J(time, x, y, z)) *
                   input(i, y_index, k)) +
                  (gamma_1 * G(time, x, y, z) * input(i + 1, y_index, k)) +
                  (P_3(time, x, y, z) * input(i, y_index, k + 1)) +
                  (four * gamma_2 * H(time, x, y, z) * input(i, y_index + 1, k)) -
                  (gamma_2 * H(time, x, y, z) * input(i, y_index + 2, k));
            solution(i, k, val);
        }
    }
}

void explicit_hhw_boundary_scheme::rhs_source(hhw_implicit_coefficients_ptr const &cfg,
                                              grid_config_3d_ptr const &grid_cfg, std::size_t const &y_index,
                                              double const &y, boundary_3d_pair const &x_boundary_pair,
                                              boundary_3d_pair const &z_boundary_pair,
                                              container_3d<by_enum::RowPlane> const &input,
                                              container_3d<by_enum::RowPlane> const &inhom_input, double const &time,
                                              container_2d<by_enum::Row> &solution)
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
        container_t lower(L, double{0.0});
        for (std::size_t k = 0; k < L + 1; ++k)
        {
            z = grid_3d::value_3(grid_cfg, k);
            lower[k] = fun(z);
        }
        solution(0, lower);
    }
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_x_bnd))
    {
        x = grid_3d::value_1(grid_cfg, N);
        auto fun = [&](double z) { return two * h_1 * ptr->value(time, y, z); };
        container_t upper(L, double{0.0});
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
        solution(N, upper);
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

hhw_explicit_boundary_solver::hhw_explicit_boundary_solver(hhw_implicit_coefficients_ptr const &coefficients,
                                                           grid_config_3d_ptr const &grid_config)
    : coefficients_{coefficients}, grid_cfg_{grid_config}
{
}

hhw_explicit_boundary_solver::~hhw_explicit_boundary_solver()
{
}

void hhw_explicit_boundary_solver::solve(container_3d<by_enum::RowPlane> const &prev_solution,
                                         boundary_3d_pair const &x_boundary_pair,
                                         boundary_3d_ptr const &y_upper_boundary_ptr,
                                         boundary_3d_pair const &z_boundary_pair, double const &time,
                                         container_3d<by_enum::RowPlane> &solution)
{
    container_3d<by_enum::ColumnPlane> cpsolution(solution);
    // 2D container for intermediate solution:
    container_2d<by_enum::Row> solution_v(coefficients_->space_size_x_, coefficients_->space_size_z_, double{});
    /// get the right-hand side of the scheme:
    auto const y = grid_3d::value_2(grid_cfg_, 0);
    explicit_hhw_boundary_scheme::rhs(coefficients_, grid_cfg_, 0, y, x_boundary_pair, z_boundary_pair, prev_solution,
                                      time, solution_v);
    // boundary plane at Y = 0:
    cpsolution(0, solution_v);
    auto const &upper_bnd_ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(y_upper_boundary_ptr);
    auto const &upper_bnd = [=](double t, double s, double r) { return upper_bnd_ptr->value(t, s, r); };
    d_2d::of_function(grid_cfg_->grid_13(), time, upper_bnd, solution_v);

    // boundary plane at Y = M-1:
    cpsolution(coefficients_->space_size_y_ - 1, solution_v);
    solution = cpsolution;
}

void hhw_explicit_boundary_solver::solve(container_3d<by_enum::RowPlane> const &prev_solution,
                                         boundary_3d_pair const &x_boundary_pair,
                                         boundary_3d_pair const &z_boundary_pair, double const &time,
                                         container_3d<by_enum::RowPlane> &solution)
{
    // 3D container for intermediate solution:
    container_2d<by_enum::Row> solution_yz(coefficients_->space_size_y_, coefficients_->space_size_z_, double{});
    container_2d<by_enum::Row> solution_xy(coefficients_->space_size_x_, coefficients_->space_size_y_, double{});
    // some constants:
    auto const &start_y = coefficients_->rangey_->lower();
    // prepare grid_xy:
    auto const &grid_12 = grid_cfg_->grid_12();
    // prepare grid_zy:
    auto const &grid_23 = grid_cfg_->grid_23();
    /// populating lower X:
    auto const &lower_x_ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(x_boundary_pair.first);
    auto const &lower_x_bnd = [=](double t, double v, double r) { return lower_x_ptr->value(t, v, r); };
    d_2d::of_function(grid_23, time, lower_x_bnd, solution_yz);
    solution(0, solution_yz);

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
    solution(lri, solution_yz);

    // populating lower Z:
    container_3d<by_enum::LayerPlane> layer_solution(solution);
    auto const &lower_z_ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(z_boundary_pair.first);
    auto const &lower_z_bnd = [=](double t, double s, double v) {
        const std::size_t i = grid_3d::index_of_1(grid_cfg_, s);
        const std::size_t j = grid_3d::index_of_2(grid_cfg_, v);
        auto const bnd_val = lower_z_ptr->value(t, s, v);
        auto const h_3 = grid_3d::step_3(grid_cfg_);
        return (((four * solution(i, j, 1)) - solution(i, j, 2) + (two * h_3 * bnd_val)) / three);
    };
    d_2d::of_function(grid_12, time, lower_z_bnd, solution_xy);
    layer_solution(0, solution_xy);

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
    layer_solution(lli, solution_xy);
    solution = layer_solution;
}

} // namespace three_dimensional

} // namespace lss_pde_solvers
