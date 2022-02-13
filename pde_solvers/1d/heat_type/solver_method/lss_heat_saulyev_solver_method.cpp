#include "lss_heat_saulyev_solver_method.hpp"

#include "../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_grids::grid_1d;
using d_1d = discretization_1d;

void heat_saulyev_solver_method::initialize(bool is_heat_sourse_set)
{
    auto a = coefficients_->A_;
    auto b = coefficients_->B_;
    auto d = coefficients_->D_;
    auto K = coefficients_->K_;

    up_sweeper_ = [=](container_t &up_component, container_t const &rhs, double time, double rhs_coeff) {
        double x{};
        for (std::size_t t = 1; t < up_component.size() - 1; ++t)
        {
            x = grid_1d::value(grid_cfg_, t);
            up_component[t] = b(time, x) * up_component[t] + d(time, x) * up_component[t + 1] +
                              a(time, x) * up_component[t - 1] + K(time, x) * rhs_coeff * rhs[t];
        }
    };

    down_sweeper_ = [=](container_t &down_component, container_t const &rhs, double time, double rhs_coeff) {
        double x{};
        for (std::size_t t = down_component.size() - 2; t >= 1; --t)
        {
            x = grid_1d::value(grid_cfg_, t);
            down_component[t] = b(time, x) * down_component[t] + d(time, x) * down_component[t + 1] +
                                a(time, x) * down_component[t - 1] + K(time, x) * rhs_coeff * rhs[t];
        }
    };

    if (is_heat_sourse_set)
    {
        source_.resize(coefficients_->space_size_);
        source_next_.resize(coefficients_->space_size_);
    }
    else
    {
        source_dummy_.resize(coefficients_->space_size_);
    }
}

heat_saulyev_solver_method::heat_saulyev_solver_method(heat_saulyev_coefficients_ptr const &coefficients,
                                                       grid_config_1d_ptr const &grid_config, bool is_heat_source_set)
    : coefficients_{coefficients}, heat_explicit_solver_method(grid_config)
{
    initialize(is_heat_source_set);
}

heat_saulyev_solver_method ::~heat_saulyev_solver_method()
{
}

void heat_saulyev_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                       std::size_t const &time_idx, double const &time, container_t &solution)
{
    if (time_idx % 2 == 0)
    {
        down_sweeper_(prev_solution, source_dummy_, time, czero_);
    }
    else
    {
        up_sweeper_(prev_solution, source_dummy_, time, czero_);
    }
    solution = prev_solution;
    solution[0] = boundary_pair.first->value(time);
    solution[coefficients_->space_size_ - 1] = boundary_pair.second->value(time);
}

void heat_saulyev_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                       std::size_t const &time_idx, double const &time, double const &next_time,
                                       std::function<double(double, double)> const &heat_source, container_t &solution)
{
    if (time_idx % 2 == 0)
    {
        d_1d::of_function(grid_cfg_, time, heat_source, source_);
        down_sweeper_(prev_solution, source_, time, cone_);
    }
    else
    {
        d_1d::of_function(grid_cfg_, next_time, heat_source, source_next_);
        up_sweeper_(prev_solution, source_next_, time, cone_);
    }
    solution = prev_solution;
    solution[0] = boundary_pair.first->value(time);
    solution[coefficients_->space_size_ - 1] = boundary_pair.second->value(time);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
