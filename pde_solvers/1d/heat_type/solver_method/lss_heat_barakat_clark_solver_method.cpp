#include "lss_heat_barakat_clark_solver_method.hpp"

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

void heat_barakat_clark_solver_method::initialize(bool is_heat_sourse_set)
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

heat_barakat_clark_solver_method::heat_barakat_clark_solver_method(
    heat_barakat_clark_coefficients_ptr const &coefficients, grid_config_1d_ptr const &grid_config,
    bool is_heat_sourse_set)
    : coefficients_{coefficients}, heat_explicit_solver_method(grid_config)
{
    initialize(is_heat_sourse_set);
}

heat_barakat_clark_solver_method::~heat_barakat_clark_solver_method()
{
}

void heat_barakat_clark_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                             std::size_t const &time_idx, double const &time, container_t &solution)
{
    // solution size:
    const std::size_t sol_size = prev_solution.size();
    //  components of the solution:
    container_t cont_up(prev_solution);
    container_t cont_down(prev_solution);
    std::future<void> sweep_up =
        std::async(std::launch::async, up_sweeper_, std::ref(cont_up), source_dummy_, time, czero_);
    std::future<void> sweep_down =
        std::async(std::launch::async, down_sweeper_, std::ref(cont_down), source_dummy_, time, czero_);
    sweep_up.wait();
    sweep_down.wait();
    cont_up[0] = cont_down[0] = boundary_pair.first->value(time);
    cont_up[sol_size - 1] = cont_down[sol_size - 1] = boundary_pair.second->value(time);
    for (std::size_t t = 0; t < sol_size; ++t)
    {
        solution[t] = chalf_ * (cont_up[t] + cont_down[t]);
    }
}

void heat_barakat_clark_solver_method::solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair,
                                             std::size_t const &time_idx, double const &time, double const &next_time,
                                             std::function<double(double, double)> const &heat_source,
                                             container_t &solution)
{
    // solution size:
    const std::size_t sol_size = prev_solution.size();
    //  components of the solution:
    container_t cont_up(prev_solution);
    container_t cont_down(prev_solution);
    d_1d::of_function(grid_cfg_, time, heat_source, source_);
    d_1d::of_function(grid_cfg_, next_time, heat_source, source_next_);
    std::future<void> sweep_up =
        std::async(std::launch::async, up_sweeper_, std::ref(cont_up), source_next_, time, cone_);
    std::future<void> sweep_down =
        std::async(std::launch::async, down_sweeper_, std::ref(cont_down), source_, time, cone_);
    sweep_up.wait();
    sweep_down.wait();
    cont_up[0] = cont_down[0] = boundary_pair.first->value(time);
    cont_up[sol_size - 1] = cont_down[sol_size - 1] = boundary_pair.second->value(time);
    for (std::size_t t = 0; t < sol_size; ++t)
    {
        solution[t] = chalf_ * (cont_up[t] + cont_down[t]);
    }
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
