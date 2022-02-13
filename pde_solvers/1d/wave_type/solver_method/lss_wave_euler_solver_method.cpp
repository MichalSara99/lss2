#include "lss_wave_euler_solver_method.hpp"

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

void explicit_wave_scheme::rhs(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                               container_t const &input_0, container_t const &input_1,
                               boundary_1d_pair const &boundary_pair, double const &time, container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
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
        solution[0] = beta * a(time, x) + c(time, x) * input_1[0] + (a(time, x) + b(time, x)) * input_1[1] -
                      d(time, x) * input_0[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * a(time, x) + (c(time, x) + alpha * a(time, x)) * input_1[0] +
                      (a(time, x) + b(time, x)) * input_1[1] - d(time, x) * input_0[0];
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
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + c(time, x) * input_1[N] - delta * b(time, x) -
                      d(time, x) * input_0[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + (c(time, x) - gamma * b(time, x)) * input_1[N] -
                      delta * b(time, x) - d(time, x) * input_0[N];
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (a(time, x) * input_1[t - 1]) + (c(time, x) * input_1[t]) + (b(time, x) * input_1[t + 1]) -
                      (d(time, x) * input_0[t]);
    }
}

void explicit_wave_scheme::rhs_source(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                                      container_t const &input_0, container_t const &input_1,
                                      container_t const &inhom_input, boundary_1d_pair const &boundary_pair,
                                      double const &time, container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
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
        solution[0] = beta * a(time, x) + c(time, x) * input_1[0] + (a(time, x) + b(time, x)) * input_1[1] -
                      d(time, x) * input_0[0] + inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * a(time, x) + (c(time, x) + alpha * a(time, x)) * input_1[0] +
                      (a(time, x) + b(time, x)) * input_1[1] - d(time, x) * input_0[0] + inhom_input[0];
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
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + c(time, x) * input_1[N] - delta * b(time, x) -
                      d(time, x) * input_0[N] + inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (a(time, x) + b(time, x)) * input_1[N - 1] + (c(time, x) - gamma * b(time, x)) * input_1[N] -
                      delta * b(time, x) - d(time, x) * input_0[N] + inhom_input[N];
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (a(time, x) * input_1[t - 1]) + (c(time, x) * input_1[t]) + (b(time, x) * input_1[t + 1]) -
                      (d(time, x) * input_0[t]) + inhom_input[t];
    }
}

void explicit_wave_scheme::rhs_initial(wave_explicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                                       container_t const &input_0, container_t const &input_1,
                                       boundary_1d_pair const &boundary_pair, double const &time, container_t &solution)
{
    const double one = 1.0;
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
    auto const &d = cfs->D_;
    auto const h = grid_1d::step(grid_config);
    auto const k = cfs->k_;
    auto const one_gamma = (two * k);
    auto const &defl = [&](double t, double x) { return (one + d(t, x)); };
    auto const &A = [&](double t, double x) { return (a(t, x) / defl(t, x)); };
    auto const &B = [&](double t, double x) { return (b(t, x) / defl(t, x)); };
    auto const &C = [&](double t, double x) { return (c(t, x) / defl(t, x)); };
    auto const &D = [&](double t, double x) { return (one_gamma * (d(t, x) / defl(t, x))); };

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
        solution[0] = beta * A(time, x) + (C(time, x) * input_0[0]) + ((A(time, x) + B(time, x)) * input_0[1]) +
                      (D(time, x) * input_1[0]);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * A(time, x) + (C(time, x) + alpha * A(time, x)) * input_0[0] +
                      (A(time, x) + B(time, x)) * input_0[1] + (D(time, x) * input_1[0]);
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
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + C(time, x) * input_0[N] - delta * B(time, x) +
                      (D(time, x) * input_1[N]);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + (C(time, x) - gamma * B(time, x)) * input_0[N] -
                      delta * B(time, x) + (D(time, x) * input_1[N]);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (A(time, x) * input_0[t - 1]) + (C(time, x) * input_0[t]) + (B(time, x) * input_0[t + 1]) +
                      (D(time, x) * input_1[t]);
    }
}

void explicit_wave_scheme::rhs_initial_source(wave_explicit_coefficients_ptr const &cfs,
                                              grid_config_1d_ptr const &grid_config, container_t const &input_0,
                                              container_t const &input_1, container_t const &inhom_input,
                                              boundary_1d_pair const &boundary_pair, double const &time,
                                              container_t &solution)
{
    const double one = 1.0;
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
    auto const &d = cfs->D_;
    auto const h = grid_1d::step(grid_config);
    auto const k = cfs->k_;
    auto const &defl = [=](double t, double x) { return (one + d(t, x)); };
    auto const &A = [=](double t, double x) { return (a(t, x) / defl(t, x)); };
    auto const &B = [=](double t, double x) { return (b(t, x) / defl(t, x)); };
    auto const &C = [=](double t, double x) { return (c(t, x) / defl(t, x)); };
    auto const one_gamma = (two * k);

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
        solution[0] = beta * A(time, x) + C(time, x) * input_0[0] + (A(time, x) + B(time, x)) * input_0[1] +
                      ((one_gamma * d(time, x)) / defl(time, x)) * input_1[0] + (inhom_input[0] / defl(time, x));
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * A(time, x) + (C(time, x) + alpha * A(time, x)) * input_0[0] +
                      (A(time, x) + B(time, x)) * input_0[1] + ((one_gamma * d(time, x)) / defl(time, x)) * input_1[0] +
                      (inhom_input[0] / defl(time, x));
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
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + C(time, x) * input_0[N] - delta * B(time, x) +
                      ((one_gamma * d(time, x)) / defl(time, x)) * input_1[N] + (inhom_input[N] / defl(time, x));
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + (C(time, x) - gamma * B(time, x)) * input_0[N] -
                      delta * B(time, x) + ((one_gamma * d(time, x)) / defl(time, x)) * input_1[N] +
                      (inhom_input[N] / defl(time, x));
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (A(time, x) * input_0[t - 1]) + (C(time, x) * input_0[t]) + (B(time, x) * input_0[t + 1]) +
                      ((one_gamma * d(time, x) / defl(time, x)) * input_1[t]) + (inhom_input[t] / defl(time, x));
    }
}

void explicit_wave_scheme::rhs_terminal(wave_explicit_coefficients_ptr const &cfs,
                                        grid_config_1d_ptr const &grid_config, container_t const &input_0,
                                        container_t const &input_1, boundary_1d_pair const &boundary_pair,
                                        double const &time, container_t &solution)
{
    const double one = 1.0;
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
    auto const &d = cfs->D_;
    auto const h = grid_1d::step(grid_config);
    auto const k = cfs->k_;
    auto const &defl = [=](double t, double x) { return (one + d(t, x)); };
    auto const &A = [=](double t, double x) { return (a(t, x) / defl(t, x)); };
    auto const &B = [=](double t, double x) { return (b(t, x) / defl(t, x)); };
    auto const &C = [=](double t, double x) { return (c(t, x) / defl(t, x)); };
    auto const one_gamma = (two * k);

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
        solution[0] = beta * A(time, x) + C(time, x) * input_0[0] + (A(time, x) + B(time, x)) * input_0[1] -
                      ((one_gamma * d(time, x)) / defl(time, x)) * input_1[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * A(time, x) + (C(time, x) + alpha * A(time, x)) * input_0[0] +
                      (A(time, x) + B(time, x)) * input_0[1] - ((one_gamma * d(time, x)) / defl(time, x)) * input_1[0];
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
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + C(time, x) * input_0[N] - delta * B(time, x) -
                      ((one_gamma * d(time, x)) / defl(time, x)) * input_1[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + (C(time, x) - gamma * B(time, x)) * input_0[N] -
                      delta * B(time, x) - ((one_gamma * d(time, x)) / defl(time, x)) * input_1[N];
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (A(time, x) * input_0[t - 1]) + (C(time, x) * input_0[t]) + (B(time, x) * input_0[t + 1]) -
                      ((one_gamma * d(time, x) / defl(time, x)) * input_1[t]);
    }
}

void explicit_wave_scheme::rhs_terminal_source(wave_explicit_coefficients_ptr const &cfs,
                                               grid_config_1d_ptr const &grid_config, container_t const &input_0,
                                               container_t const &input_1, container_t const &inhom_input,
                                               boundary_1d_pair const &boundary_pair, double const &time,
                                               container_t &solution)
{
    const double one = 1.0;
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &a = cfs->A_;
    auto const &b = cfs->B_;
    auto const &c = cfs->C_;
    auto const &d = cfs->D_;
    auto const h = grid_1d::step(grid_config);
    auto const k = cfs->k_;
    auto const &defl = [=](double t, double x) { return (one + d(t, x)); };
    auto const &A = [=](double t, double x) { return (a(t, x) / defl(t, x)); };
    auto const &B = [=](double t, double x) { return (b(t, x) / defl(t, x)); };
    auto const &C = [=](double t, double x) { return (c(t, x) / defl(t, x)); };
    auto const one_gamma = (two * k);

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
        solution[0] = beta * A(time, x) + C(time, x) * input_0[0] + (A(time, x) + B(time, x)) * input_0[1] -
                      ((one_gamma * d(time, x)) / defl(time, x)) * input_1[0] + (inhom_input[0] / defl(time, x));
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = beta * A(time, x) + (C(time, x) + alpha * A(time, x)) * input_0[0] +
                      (A(time, x) + B(time, x)) * input_0[1] - ((one_gamma * d(time, x)) / defl(time, x)) * input_1[0] +
                      (inhom_input[0] / defl(time, x));
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
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + C(time, x) * input_0[N] - delta * B(time, x) -
                      ((one_gamma * d(time, x)) / defl(time, x)) * input_1[N] + (inhom_input[N] / defl(time, x));
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (A(time, x) + B(time, x)) * input_0[N - 1] + (C(time, x) - gamma * B(time, x)) * input_0[N] -
                      delta * B(time, x) - ((one_gamma * d(time, x)) / defl(time, x)) * input_1[N] +
                      (inhom_input[N] / defl(time, x));
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_config, t);
        solution[t] = (A(time, x) * input_0[t - 1]) + (C(time, x) * input_0[t]) + (B(time, x) * input_0[t + 1]) -
                      ((one_gamma * d(time, x) / defl(time, x)) * input_1[t]) + (inhom_input[t] / defl(time, x));
    }
}

void wave_euler_solver_method::initialize(bool is_wave_source_set)
{
    if (is_wave_source_set)
    {
        source_.resize(coefficients_->space_size_);
    }
}

wave_euler_solver_method::wave_euler_solver_method(wave_explicit_coefficients_ptr const &coefficients,
                                                   grid_config_1d_ptr const &grid_config, bool is_wave_source_set)
    : coefficients_{coefficients}, wave_explicit_solver_method(grid_config)
{
    initialize(is_wave_source_set);
}

wave_euler_solver_method ::~wave_euler_solver_method()
{
}

void wave_euler_solver_method::solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                                             boundary_1d_pair const &boundary_pair, double const &time,
                                             double const &next_time, container_t &solution)
{
    // get the right-hand side of the scheme:
    explicit_wave_scheme::rhs_initial(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                                      solution);
}

void wave_euler_solver_method::solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                                             boundary_1d_pair const &boundary_pair, double const &time,
                                             double const &next_time,
                                             std::function<double(double, double)> const &wave_source,
                                             container_t &solution)
{
    // get the right-hand side of the scheme:
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    explicit_wave_scheme::rhs_initial_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                             boundary_pair, time, solution);
}

void wave_euler_solver_method::solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                              boundary_1d_pair const &boundary_pair, double const &time,
                                              double const &next_time, container_t &solution)
{
    // get the right-hand side of the scheme:
    explicit_wave_scheme::rhs_terminal(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                                       solution);
}

void wave_euler_solver_method::solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                              boundary_1d_pair const &boundary_pair, double const &time,
                                              double const &next_time,
                                              std::function<double(double, double)> const &wave_source,
                                              container_t &solution)
{
    // get the right-hand side of the scheme:
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    explicit_wave_scheme::rhs_terminal_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                              boundary_pair, time, solution);
}

void wave_euler_solver_method::solve(container_t &prev_solution_0, container_t &prev_solution_1,
                                     boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                                     container_t &solution)
{
    // get the right-hand side of the scheme:
    explicit_wave_scheme::rhs(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                              solution);
}

void wave_euler_solver_method::solve(container_t &prev_solution_0, container_t &prev_solution_1,
                                     boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                                     std::function<double(double, double)> const &wave_source, container_t &solution)
{
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    explicit_wave_scheme::rhs_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_, boundary_pair,
                                     time, solution);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
