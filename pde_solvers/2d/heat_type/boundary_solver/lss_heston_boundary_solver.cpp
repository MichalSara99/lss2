#include "lss_heston_boundary_solver.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../lss_pde_discretization_config.hpp"

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_2d;
using lss_boundary::neumann_boundary_2d;
using d_1d = discretization_1d;
using lss_grids::grid_2d;

namespace two_dimensional
{

void heston_boundary_rhs::rhs(heat_coefficients_2d_ptr const &cfg, grid_config_2d_ptr const &grid_cfg,
                              std::size_t const &y_index, double const &y,
                              boundary_2d_pair const &horizontal_boundary_pair, matrix_2d const &input,
                              double const &time, container_t &solution)
{
    auto const four = 4.0;
    auto const three = 3.0;
    auto const two = 2.0;
    auto const one = 1.0;
    auto const &D = cfg->D_;
    auto const &E = cfg->E_;
    auto const &F = cfg->F_;
    auto const delta = cfg->delta_;
    auto const ni = cfg->ni_;
    auto const rho = cfg->rho_;

    auto const &first_bnd = horizontal_boundary_pair.first;
    auto const &second_bnd = horizontal_boundary_pair.second;

    double x{}, h_1{};
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(first_bnd))
    {
        solution[0] = ptr->value(time, y);
    }

    const std::size_t N = solution.size() - 1;
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(second_bnd))
    {
        x = grid_2d::value_1(grid_cfg, N);
        h_1 = grid_2d::step_1(grid_cfg);
        const double delta_ = two * h_1 * ptr->value(time, y);
        solution[N] = ((one - three * ni * E(time, x, y) + rho * F(time, x, y)) * input(N, y_index)) +
                      (four * ni * E(time, x, y) * input(N, y_index + 1)) -
                      (ni * E(time, x, y) * input(N, y_index + 2)) - (delta * delta_ * D(time, x, y));
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        solution[N] = ptr->value(time, y);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_2d::value_1(grid_cfg, t);
        solution[t] = (-delta * D(time, x, y) * input(t - 1, y_index)) +
                      ((one - three * ni * E(time, x, y) + rho * F(time, x, y)) * input(t, y_index)) +
                      (delta * D(time, x, y) * input(t + 1, y_index)) - (ni * E(time, x, y) * input(t, y_index + 2)) +
                      (four * ni * E(time, x, y) * input(t, y_index + 1));
    }
}

void heston_boundary_rhs::rhs_source(heat_coefficients_2d_ptr const &cfg, grid_config_2d_ptr const &grid_cfg,
                                     std::size_t const &y_index, double const &y,
                                     boundary_2d_pair const &horizontal_boundary_pair, matrix_2d const &input,
                                     matrix_2d const &inhom_input, double const &time, container_t &solution)
{
    auto const four = 4.0;
    auto const three = 3.0;
    auto const two = 2.0;
    auto const one = 1.0;
    auto const &D = cfg->D_;
    auto const &E = cfg->E_;
    auto const &F = cfg->F_;
    auto const delta = cfg->delta_;
    auto const ni = cfg->ni_;
    auto const rho = cfg->rho_;

    auto const &first_bnd = horizontal_boundary_pair.first;
    auto const &second_bnd = horizontal_boundary_pair.second;

    double x{}, h_1{};
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(first_bnd))
    {
        solution[0] = ptr->value(time, y);
    }

    const std::size_t N = solution.size() - 1;
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(second_bnd))
    {
        x = grid_2d::value_1(grid_cfg, N);
        h_1 = grid_2d::step_1(grid_cfg);
        const double delta_ = two * h_1 * ptr->value(time, y);
        solution[N] = ((one - three * ni * E(time, x, y) + rho * F(time, x, y)) * input(N, y_index)) +
                      (four * ni * E(time, x, y) * input(N, y_index + 1)) -
                      (ni * E(time, x, y) * input(N, y_index + 2)) - (delta * delta_ * D(time, x, y)) +
                      (rho * inhom_input(N, y_index));
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        solution[N] = ptr->value(time, y);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_2d::value_1(grid_cfg, t);
        solution[t] = (-delta * D(time, x, y) * input(t - 1, y_index)) +
                      ((one - three * ni * E(time, x, y) + rho * F(time, x, y)) * input(t, y_index)) +
                      (delta * D(time, x, y) * input(t + 1, y_index)) - (ni * E(time, x, y) * input(t, y_index + 2)) +
                      (four * ni * E(time, x, y) * input(t, y_index + 1)) + (rho * inhom_input(N, y_index));
    }
}

heston_boundary_solver::heston_boundary_solver(heat_coefficients_2d_ptr const &coefficients,
                                               grid_config_2d_ptr const &grid_config)
    : coefficients_{coefficients}, grid_cfg_{grid_config}
{
}

heston_boundary_solver::~heston_boundary_solver()
{
}

void heston_boundary_solver::solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizonatal_boundary_pair,
                                   boundary_2d_ptr const &vertical_upper_boundary_ptr, double const &time,
                                   matrix_2d &solution)
{
    // 1D container for intermediate solution:
    container_t solution_v(double{}, coefficients_->space_size_x_);
    // get the right-hand side of the scheme:
    auto const y = grid_2d::value_2(grid_cfg_, 0);
    heston_boundary_rhs::rhs(coefficients_, grid_cfg_, 0, y, horizonatal_boundary_pair, prev_solution, time,
                             solution_v);
    solution.column(0) = solution_v;
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(vertical_upper_boundary_ptr))
    {
        auto const &upper_bnd = [=](double t, double s) { return ptr->value(t, s); };
        d_1d::of_function(grid_cfg_->grid_1(), time, upper_bnd, solution_v);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(vertical_upper_boundary_ptr))
    {
        auto const &grid_1 = grid_cfg_->grid_1();
        auto const &grid_2 = grid_cfg_->grid_2();
        auto const lci = solution.columns() - 1;
        auto const &upper_bnd = [=](double t, double s) {
            const std::size_t i = grid_1d::index_of(grid_1, s);
            auto const bnd_val = ptr->value(t, s);
            auto const h_2 = grid_1d::step(grid_2);
            return (((4.0 * prev_solution(i, lci - 1)) - prev_solution(i, lci - 2) - (2.0 * h_2 * bnd_val)) / 3.0);
        };
        d_1d::of_function(grid_cfg_->grid_1(), time, upper_bnd, solution_v);
    }
    solution.column(coefficients_->space_size_y_ - 1) = solution_v;
}

void heston_boundary_solver::solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizonatal_boundary_pair,
                                   double const &time, matrix_2d &solution)
{
    // 1D container for intermediate solution:
    container_t solution_v(double{}, coefficients_->space_size_y_);
    // some constants:
    auto const &start_y = coefficients_->rangey_->lower();
    // prepare grid_1:
    auto const &grid_1 = grid_cfg_->grid_1();
    // prepare grid_2:
    auto const &grid_2 = grid_cfg_->grid_2();
    // populating lower horizontal:
    auto const &lower_bnd_ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(horizonatal_boundary_pair.first);
    auto const &lower_bnd = [=](double t, double v) { return lower_bnd_ptr->value(t, v); };
    d_1d::of_function(grid_2, time, lower_bnd, solution_v);
    solution.row(0) = solution_v;
    // populating upper horizontal:
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(horizonatal_boundary_pair.second))
    {
        auto const lri = solution.rows() - 1;
        auto const &upper_bnd = [=](double t, double v) {
            const std::size_t j = grid_1d::index_of(grid_2, v);
            auto const bnd_val = ptr->value(t, v);
            auto const h_1 = grid_1d::step(grid_1);
            return (((4.0 * solution(lri - 1, j)) - solution(lri - 2, j) - (2.0 * h_1 * bnd_val)) / 3.0);
        };
        d_1d::of_function(grid_2, time, upper_bnd, solution_v);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(horizonatal_boundary_pair.second))
    {
        auto const &upper_bnd = [=](double t, double v) { return ptr->value(t, v); };
        d_1d::of_function(grid_2, time, upper_bnd, solution_v);
    }
    auto const N_x = coefficients_->space_size_x_;
    solution.row(N_x - 1) = solution_v;
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
