#include "lss_heat_implicit_solver_method.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"
#include "../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using lss_grids::grid_1d;
using d_1d = discretization_1d;

void implicit_heat_scheme::rhs(heat_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                               container_t const &input, boundary_1d_pair const &boundary_pair, double const &time,
                               container_t &solution)
{
    auto const two = 2.0;
    auto const one = 1.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &D = cfs->D_;
    auto const theta = cfs->theta_;
    auto const h = grid_1d::step(grid_cfg);
    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = (one - theta) * beta * A(time, x) + (one - two * (one - theta) * B(time, x)) * input[0] +
                      (one - theta) * (A(time, x) + D(time, x)) * input[1];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (one - theta) * beta * A(time, x) +
                      (one - (one - theta) * (two * B(time, x) - alpha * A(time, x))) * input[0] +
                      (one - theta) * (A(time, x) + D(time, x)) * input[1];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (one - theta) * (A(time, x) + D(time, x)) * input[N - 1] +
                      (one - two * (one - theta) * B(time, x)) * input[N] - (one - theta) * delta * D(time, x);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (one - theta) * (A(time, x) + D(time, x)) * input[N - 1] +
                      (one - (one - theta) * (two * B(time, x) + gamma * D(time, x))) * input[N] -
                      (one - theta) * delta * D(time, x);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (D(time, x) * (one - theta) * input[t + 1]) +
                      ((one - two * B(time, x) * (one - theta)) * input[t]) +
                      (A(time, x) * (one - theta) * input[t - 1]);
    }
}

void implicit_heat_scheme::rhs_source(heat_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                                      container_t const &input, container_t const &inhom_input,
                                      container_t const &inhom_input_next, boundary_1d_pair const &boundary_pair,
                                      double const &time, container_t &solution)
{
    auto const two = 2.0;
    auto const one = 1.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const theta = cfs->theta_;
    auto const h = grid_1d::step(grid_cfg);
    double x{};

    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = (one - theta) * beta * A(time, x) + (one - two * (one - theta) * B(time, x)) * input[0] +
                      (one - theta) * (A(time, x) + D(time, x)) * input[1] + theta * k * inhom_input_next[0] +
                      (one - theta) * k * inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (one - theta) * beta * A(time, x) +
                      (one - (one - theta) * (two * B(time, x) - alpha * A(time, x))) * input[0] +
                      (one - theta) * (A(time, x) + D(time, x)) * input[1] + theta * k * inhom_input_next[0] +
                      (one - theta) * k * inhom_input[0];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (one - theta) * (A(time, x) + D(time, x)) * input[N - 1] +
                      (one - two * (one - theta) * B(time, x)) * input[N] - (one - theta) * delta * D(time, x) +
                      theta * k * inhom_input_next[N] + (one - theta) * k * inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (one - theta) * (A(time, x) + D(time, x)) * input[N - 1] +
                      (one - (one - theta) * (two * B(time, x) + gamma * D(time, x))) * input[N] -
                      (one - theta) * delta * D(time, x) + theta * k * inhom_input_next[N] +
                      (one - theta) * k * inhom_input[N];
    }
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (D(time, x) * (one - theta) * input[t + 1]) +
                      ((one - two * B(time, x) * (one - theta)) * input[t]) +
                      (A(time, x) * (one - theta) * input[t - 1]) +
                      k * (theta * inhom_input_next[t] + (one - theta) * inhom_input[t]);
    }
}

void heat_implicit_solver_method::initialize(bool is_heat_sourse_set)
{
    // prepare containers:
    low_.resize(coefficients_->space_size_);
    diag_.resize(coefficients_->space_size_);
    high_.resize(coefficients_->space_size_);
    rhs_.resize(coefficients_->space_size_);
    if (is_heat_sourse_set)
    {
        source_.resize(coefficients_->space_size_);
        source_next_.resize(coefficients_->space_size_);
    }
}

void heat_implicit_solver_method::split(double const &time, container_t &low, container_t &diag, container_t &high)
{
    double x{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        x = grid_1d::value(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->A_(time, x));
        diag[t] = (cone_ + ctwo_ * coefficients_->theta_ * coefficients_->B_(time, x));
        high[t] = (-coefficients_->theta_ * coefficients_->D_(time, x));
    }
}

heat_implicit_solver_method::heat_implicit_solver_method(
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solver_ptr, heat_coefficients_ptr const &coefficients,
    grid_config_1d_ptr const &grid_config, bool is_heat_sourse_set)
    : solveru_ptr_{solver_ptr}, coefficients_{coefficients}, grid_cfg_{grid_config}
{
    initialize(is_heat_sourse_set);
}

heat_implicit_solver_method::~heat_implicit_solver_method()
{
}

void heat_implicit_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                        double const &time, container_t &solution)
{
    implicit_heat_scheme::rhs(coefficients_, grid_cfg_, prev_solution, boundary_pair, time, rhs_);
    split(time, low_, diag_, high_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, time);
}

void heat_implicit_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                        double const &time, double const &next_time,
                                        std::function<double(double, double)> const &heat_source, container_t &solution)
{
    split(time, low_, diag_, high_);
    d_1d::of_function(grid_cfg_, time, heat_source, source_);
    d_1d::of_function(grid_cfg_, next_time, heat_source, source_next_);
    implicit_heat_scheme::rhs_source(coefficients_, grid_cfg_, prev_solution, source_, source_next_, boundary_pair,
                                     time, rhs_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, time);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
