#include "lss_heat_euler_solver_method.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;
using d_1d = discretization_1d;

void explicit_heat_scheme::rhs(heat_euler_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                               container_t const &input, boundary_1d_pair const &boundary_pair, double const &time,
                               container_t &solution)
{
    auto const two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &d = cfs->D_;
    auto const h = grid_1d::step(grid_config);
    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_config, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(first_bnd))
    {
        solution[0] = ptr->value(time);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = beta * a(time, x) + b(time, x) * input[0] + (a(time, x) + d(time, x)) * input[1];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] =
            (b(time, x) + alpha * a(time, x)) * input[0] + (a(time, x) + d(time, x)) * input[1] + beta * a(time, x);
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_config, N);
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        solution[N] = ptr->value(time);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (a(time, x) + d(time, x)) * input[N - 1] + b(time, x) * input[N] - delta * d(time, x);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + d(time, x)) * input[N - 1] + (b(time, x) - gamma * d(time, x)) * input[N] -
                      delta * d(time, x);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (d(time, x) * input[t + 1]) + (b(time, x) * input[t]) + (a(time, x) * input[t - 1]);
    }
}

void explicit_heat_scheme::rhs_source(heat_euler_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                                      container_t const &input, container_t const &inhom_input,
                                      boundary_1d_pair const &boundary_pair, double const &time, container_t &solution)
{
    auto const two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &d = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_config);
    double x{};

    // for lower boundaries first:
    x = grid_1d::value(grid_config, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] =
            beta * a(time, x) + b(time, x) * input[0] + (a(time, x) + d(time, x)) * input[1] + k * inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (b(time, x) + alpha * a(time, x)) * input[0] + (a(time, x) + d(time, x)) * input[1] +
                      beta * a(time, x) + k * inhom_input[0];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_config, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] =
            (a(time, x) + d(time, x)) * input[N - 1] + b(time, x) * input[N] - delta * d(time, x) + k * inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + d(time, x)) * input[N - 1] + (b(time, x) - gamma * d(time, x)) * input[N] -
                      delta * d(time, x) + k * inhom_input[N];
        ;
    }
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] =
            (d(time, x) * input[t + 1]) + (b(time, x) * input[t]) + (a(time, x) * input[t - 1]) + (k * inhom_input[t]);
    }
}

void heat_euler_solver_method::initialize(bool is_heat_source_set)
{
    if (is_heat_source_set)
    {
        source_.resize(coefficients_->space_size_);
    }
}

heat_euler_solver_method::heat_euler_solver_method(heat_euler_coefficients_ptr const &coefficients,
                                                   grid_config_1d_ptr const &grid_config, bool is_heat_source_set)
    : coefficients_{coefficients}, heat_explicit_solver_method(grid_config)
{
    initialize(is_heat_source_set);
}

heat_euler_solver_method::~heat_euler_solver_method()
{
}

void heat_euler_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                     std::size_t const &time_idx, double const &time, container_t &solution)
{
    explicit_heat_scheme::rhs(coefficients_, grid_cfg_, prev_solution, boundary_pair, time, solution);
}

void heat_euler_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                     std::size_t const &time_idx, double const &time, double const &next_time,
                                     std::function<double(double, double)> const &heat_source, container_t &solution)
{
    d_1d::of_function(grid_cfg_, time, heat_source, source_);
    explicit_heat_scheme::rhs_source(coefficients_, grid_cfg_, prev_solution, source_, boundary_pair, time, solution);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
