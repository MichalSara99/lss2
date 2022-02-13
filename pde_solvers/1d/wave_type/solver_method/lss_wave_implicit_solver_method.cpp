#include "lss_wave_implicit_solver_method.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../boundaries/lss_robin_boundary.hpp"
#include "../../../../common/lss_macros.hpp"
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

void implicit_wave_scheme::rhs(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                               container_t const &input_0, container_t const &input_1,
                               boundary_1d_pair const &boundary_pair, double const &time, container_t &solution)
{
    const double two = 2.0;
    const double one = 1.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &C = cfs->C_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_cfg);
    auto const lambda = one / (k * k);

    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta_curr = two * h * ptr->value(time);
        const double beta_prev = two * h * ptr->value(time - k);
        solution[0] = two * (A(time, x) + B(time, x)) * input_1[1] + two * (lambda - C(time, x)) * input_1[0] +
                      (A(time, x) + B(time, x)) * input_0[1] - (D(time, x) + C(time, x)) * input_0[0] +
                      (two * beta_curr + beta_prev) * A(time, x);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta_curr = two * h * ptr->value(time);
        const double beta_prev = two * h * ptr->value(time - k);
        const double alpha_curr = two * h * ptr->linear_value(time);
        const double alpha_prev = two * h * ptr->linear_value(time - k);
        solution[0] = two * (A(time, x) + B(time, x)) * input_1[1] +
                      two * (lambda - C(time, x) + alpha_curr * A(time, x)) * input_1[0] +
                      (A(time, x) + B(time, x)) * input_0[1] -
                      (D(time, x) + C(time, x) - alpha_prev * A(time, x)) * input_0[0] +
                      (two * beta_curr + beta_prev) * A(time, x);
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta_curr = two * h * ptr->value(time);
        const double delta_prev = two * h * ptr->value(time - k);
        solution[N] = two * (A(time, x) + B(time, x)) * input_1[N - 1] + two * (lambda - C(time, x)) * input_1[N] +
                      (A(time, x) + B(time, x)) * input_0[N - 1] - (D(time, x) + C(time, x)) * input_0[N] -
                      (two * delta_curr + delta_prev) * B(time, x);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta_curr = two * h * ptr->value(time);
        const double delta_prev = two * h * ptr->value(time - k);
        const double gamma_curr = two * h * ptr->linear_value(time);
        const double gamma_prev = two * h * ptr->linear_value(time - k);
        solution[N] = two * (A(time, x) + B(time, x)) * input_1[N - 1] +
                      two * (lambda - C(time, x) - gamma_curr * B(time, x)) * input_1[N] +
                      (A(time, x) + B(time, x)) * input_0[N - 1] -
                      (D(time, x) + C(time, x) + gamma_prev * B(time, x)) * input_0[N] -
                      (two * delta_curr + delta_prev) * B(time, x);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (B(time, x) * input_0[t + 1]) - ((D(time, x) + C(time, x)) * input_0[t]) +
                      (A(time, x) * input_0[t - 1]) + (two * B(time, x) * input_1[t + 1]) +
                      (two * (lambda - C(time, x)) * input_1[t]) + (two * A(time, x) * input_1[t - 1]);
    }
}

void implicit_wave_scheme::rhs_source(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                                      container_t const &input_0, container_t const &input_1,
                                      container_t const &inhom_input, boundary_1d_pair const &boundary_pair,
                                      double const &time, container_t &solution)
{
    const double two = 2.0;
    const double one = 1.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &C = cfs->C_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_cfg);
    auto const lambda = one / (k * k);
    double x{};

    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta_curr = two * h * ptr->value(time);
        const double beta_prev = two * h * ptr->value(time - k);
        solution[0] = two * (A(time, x) + B(time, x)) * input_1[1] + two * (lambda - C(time, x)) * input_1[0] +
                      (A(time, x) + B(time, x)) * input_0[1] - (D(time, x) + C(time, x)) * input_0[0] +
                      (two * beta_curr + beta_prev) * A(time, x) + inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta_curr = two * h * ptr->value(time);
        const double beta_prev = two * h * ptr->value(time - k);
        const double alpha_curr = two * h * ptr->linear_value(time);
        const double alpha_prev = two * h * ptr->linear_value(time - k);
        solution[0] = two * (A(time, x) + B(time, x)) * input_1[1] +
                      two * (lambda - C(time, x) + alpha_curr * A(time, x)) * input_1[0] +
                      (A(time, x) + B(time, x)) * input_0[1] -
                      (D(time, x) + C(time, x) - alpha_prev * A(time, x)) * input_0[0] +
                      (two * beta_curr + beta_prev) * A(time, x) + inhom_input[0];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta_curr = two * h * ptr->value(time);
        const double delta_prev = two * h * ptr->value(time - k);
        solution[N] = two * (A(time, x) + B(time, x)) * input_1[N - 1] + two * (lambda - C(time, x)) * input_1[N] +
                      (A(time, x) + B(time, x)) * input_0[N - 1] - (D(time, x) + C(time, x)) * input_0[N] -
                      (two * delta_curr + delta_prev) * B(time, x) + inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta_curr = two * h * ptr->value(time);
        const double delta_prev = two * h * ptr->value(time - k);
        const double gamma_curr = two * h * ptr->linear_value(time);
        const double gamma_prev = two * h * ptr->linear_value(time - k);
        solution[N] = two * (A(time, x) + B(time, x)) * input_1[N - 1] +
                      two * (lambda - C(time, x) - gamma_curr * B(time, x)) * input_1[N] +
                      (A(time, x) + B(time, x)) * input_0[N - 1] -
                      (D(time, x) + C(time, x) + gamma_prev * B(time, x)) * input_0[N] -
                      (two * delta_curr + delta_prev) * B(time, x) + inhom_input[N];
    }
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (B(time, x) * input_0[t + 1]) - ((D(time, x) + C(time, x)) * input_0[t]) +
                      (A(time, x) * input_0[t - 1]) + (two * B(time, x) * input_1[t + 1]) +
                      (two * (lambda - C(time, x)) * input_1[t]) + (two * A(time, x) * input_1[t - 1]) + inhom_input[t];
    }
}

void implicit_wave_scheme::rhs_initial(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                                       container_t const &input_0, container_t const &input_1,
                                       boundary_1d_pair const &boundary_pair, double const &time, container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &C = cfs->C_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_cfg);
    auto const lambda = cfs->lambda_;
    auto const one_gamma = (two * k);

    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) + (two * (lambda - C(time, x)) * input_0[0]) +
                      (two * beta * A(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[0]);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) +
                      (two * (lambda - C(time, x) + alpha * A(time, x)) * input_0[0]) + (two * beta * A(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[0]);
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) + (two * (lambda - C(time, x)) * input_0[N]) -
                      (two * delta * B(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[N]);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) +
                      (two * (lambda - C(time, x) - gamma * B(time, x)) * input_0[N]) - (two * delta * B(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[N]);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (two * B(time, x) * input_0[t + 1]) + (two * (lambda - C(time, x)) * input_0[t]) +
                      (two * A(time, x) * input_0[t - 1]) - (one_gamma * B(time, x) * input_1[t + 1]) +
                      (one_gamma * (D(time, x) + C(time, x)) * input_1[t]) - (one_gamma * A(time, x) * input_1[t - 1]);
    }
}

void implicit_wave_scheme::rhs_initial_source(wave_implicit_coefficients_ptr const &cfs,
                                              grid_config_1d_ptr const &grid_cfg, container_t const &input_0,
                                              container_t const &input_1, container_t const &inhom_input,
                                              boundary_1d_pair const &boundary_pair, double const &time,
                                              container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &C = cfs->C_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_cfg);
    auto const lambda = cfs->lambda_;
    auto const one_gamma = (two * k);

    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) + (two * (lambda - C(time, x)) * input_0[0]) +
                      (two * beta * A(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[0]) + inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) +
                      (two * (lambda - C(time, x) + alpha * A(time, x)) * input_0[0]) + (two * beta * A(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[0]) + inhom_input[0];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) + (two * (lambda - C(time, x)) * input_0[N]) -
                      (two * delta * B(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[N]) + inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) +
                      (two * (lambda - C(time, x) - gamma * B(time, x)) * input_0[N]) - (two * delta * B(time, x)) +
                      (one_gamma * (D(time, x) + C(time, x) - A(time, x) - B(time, x)) * input_1[N]) + inhom_input[N];
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (two * B(time, x) * input_0[t + 1]) + (two * (lambda - C(time, x)) * input_0[t]) +
                      (two * A(time, x) * input_0[t - 1]) - (one_gamma * B(time, x) * input_1[t + 1]) +
                      (one_gamma * (D(time, x) + C(time, x)) * input_1[t]) - (one_gamma * A(time, x) * input_1[t - 1]) +
                      inhom_input[t];
    }
}

void implicit_wave_scheme::rhs_terminal(wave_implicit_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_cfg,
                                        container_t const &input_0, container_t const &input_1,
                                        boundary_1d_pair const &boundary_pair, double const &time,
                                        container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &C = cfs->C_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_cfg);
    auto const lambda = cfs->lambda_;
    auto const one_gamma = (two * k);

    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) + (two * (lambda - C(time, x)) * input_0[0]) +
                      (two * beta * A(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[0]);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) +
                      (two * (lambda - C(time, x) + alpha * A(time, x)) * input_0[0]) + (two * beta * A(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[0]);
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) + (two * (lambda - C(time, x)) * input_0[N]) -
                      (two * delta * B(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[N]);
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) +
                      (two * (lambda - C(time, x) - gamma * B(time, x)) * input_0[N]) - (two * delta * B(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[N]);
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (two * B(time, x) * input_0[t + 1]) + (two * (lambda - C(time, x)) * input_0[t]) +
                      (two * A(time, x) * input_0[t - 1]) + (one_gamma * B(time, x) * input_1[t + 1]) -
                      (one_gamma * (D(time, x) + C(time, x)) * input_1[t]) + (one_gamma * A(time, x) * input_1[t - 1]);
    }
}

void implicit_wave_scheme::rhs_terminal_source(wave_implicit_coefficients_ptr const &cfs,
                                               grid_config_1d_ptr const &grid_cfg, container_t const &input_0,
                                               container_t const &input_1, container_t const &inhom_input,
                                               boundary_1d_pair const &boundary_pair, double const &time,
                                               container_t &solution)
{
    const double two = 2.0;
    auto const &first_bnd = boundary_pair.first;
    auto const &second_bnd = boundary_pair.second;
    auto const &A = cfs->A_;
    auto const &B = cfs->B_;
    auto const &C = cfs->C_;
    auto const &D = cfs->D_;
    auto const k = cfs->k_;
    auto const h = grid_1d::step(grid_cfg);
    auto const lambda = cfs->lambda_;
    auto const one_gamma = (two * k);

    double x{};
    // for lower boundaries first:
    x = grid_1d::value(grid_cfg, 0);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) + (two * (lambda - C(time, x)) * input_0[0]) +
                      (two * beta * A(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[0]) + inhom_input[0];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const double beta = two * h * ptr->value(time);
        const double alpha = two * h * ptr->linear_value(time);
        solution[0] = (two * (A(time, x) + B(time, x)) * input_0[1]) +
                      (two * (lambda - C(time, x) + alpha * A(time, x)) * input_0[0]) + (two * beta * A(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[0]) + inhom_input[0];
    }
    // for upper boundaries second:
    const std::size_t N = solution.size() - 1;
    x = grid_1d::value(grid_cfg, N);
    if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) + (two * (lambda - C(time, x)) * input_0[N]) -
                      (two * delta * B(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[N]) + inhom_input[N];
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const double delta = two * h * ptr->value(time);
        const double gamma = two * h * ptr->linear_value(time);
        solution[N] = (two * (A(time, x) + B(time, x)) * input_0[N - 1]) +
                      (two * (lambda - C(time, x) - gamma * B(time, x)) * input_0[N]) - (two * delta * B(time, x)) +
                      (one_gamma * (A(time, x) + B(time, x) - C(time, x) - D(time, x)) * input_1[N]) + inhom_input[N];
    }

    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_1d::value(grid_cfg, t);
        solution[t] = (two * B(time, x) * input_0[t + 1]) + (two * (lambda - C(time, x)) * input_0[t]) +
                      (two * A(time, x) * input_0[t - 1]) + (one_gamma * B(time, x) * input_1[t + 1]) -
                      (one_gamma * (D(time, x) + C(time, x)) * input_1[t]) + (one_gamma * A(time, x) * input_1[t - 1]) +
                      inhom_input[t];
    }
}

void wave_implicit_solver_method::initialize(bool is_wave_source_set)
{
    low_.resize(coefficients_->space_size_);
    diag_.resize(coefficients_->space_size_);
    high_.resize(coefficients_->space_size_);
    rhs_.resize(coefficients_->space_size_);
    if (is_wave_source_set)
    {
        source_.resize(coefficients_->space_size_);
    }
}

void wave_implicit_solver_method::split_0(double time, container_t &low_0, container_t &diag_0, container_t &high_0)
{
    double x{};
    for (std::size_t t = 0; t < low_0.size(); ++t)
    {
        x = grid_1d::value(grid_cfg_, t);
        low_0[t] = ctwo_ * (-cone_ * coefficients_->A_(time, x));
        diag_0[t] = (coefficients_->E_(time, x) + ctwo_ * coefficients_->C_(time, x) + coefficients_->D_(time, x));
        high_0[t] = ctwo_ * (-cone_ * coefficients_->B_(time, x));
    }
}

void wave_implicit_solver_method::split_1(double time, container_t &low_1, container_t &diag_1, container_t &high_1)
{
    double x{};
    for (std::size_t t = 0; t < low_1.size(); ++t)
    {
        x = grid_1d::value(grid_cfg_, t);
        low_1[t] = -cone_ * coefficients_->A_(time, x);
        diag_1[t] = (coefficients_->E_(time, x) + coefficients_->C_(time, x));
        high_1[t] = -cone_ * coefficients_->B_(time, x);
    }
}

wave_implicit_solver_method ::wave_implicit_solver_method(
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solver_ptr,
    wave_implicit_coefficients_ptr const &coefficients, grid_config_1d_ptr const &grid_config, bool is_wave_source_set)
    : solveru_ptr_{solver_ptr}, coefficients_{coefficients}, grid_cfg_{grid_config}
{
    initialize(is_wave_source_set);
}

wave_implicit_solver_method::~wave_implicit_solver_method()
{
}

void wave_implicit_solver_method::solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                                                boundary_1d_pair const &boundary_pair, double const &time,
                                                double const &next_time, container_t &solution)
{
    implicit_wave_scheme::rhs_initial(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                                      rhs_);
    split_0(time, low_, diag_, high_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, next_time);
}

void wave_implicit_solver_method::solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                                                boundary_1d_pair const &boundary_pair, double const &time,
                                                double const &next_time,
                                                std::function<double(double, double)> const &wave_source,
                                                container_t &solution)
{
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    implicit_wave_scheme::rhs_initial_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                             boundary_pair, time, rhs_);
    split_0(time, low_, diag_, high_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, next_time);
}

void wave_implicit_solver_method::solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                                 boundary_1d_pair const &boundary_pair, double const &time,
                                                 double const &next_time, container_t &solution)
{
    implicit_wave_scheme::rhs_terminal(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time,
                                       rhs_);
    split_0(time, low_, diag_, high_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, next_time);
}

void wave_implicit_solver_method::solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                                 boundary_1d_pair const &boundary_pair, double const &time,
                                                 double const &next_time,
                                                 std::function<double(double, double)> const &wave_source,
                                                 container_t &solution)
{
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    implicit_wave_scheme::rhs_terminal_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_,
                                              boundary_pair, time, rhs_);
    split_0(time, low_, diag_, high_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, next_time);
}

void wave_implicit_solver_method::solve(container_t &prev_solution_0, container_t &prev_solution_1,
                                        boundary_1d_pair const &boundary_pair, double const &time,
                                        double const &next_time, container_t &solution)
{
    split_1(time, low_, diag_, high_);
    implicit_wave_scheme::rhs(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, boundary_pair, time, rhs_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, next_time);
}

void wave_implicit_solver_method::solve(container_t &prev_solution_0, container_t &prev_solution_1,
                                        boundary_1d_pair const &boundary_pair, double const &time,
                                        double const &next_time,
                                        std::function<double(double, double)> const &wave_source, container_t &solution)
{
    split_1(time, low_, diag_, high_);
    d_1d::of_function(grid_cfg_, time, wave_source, source_);
    implicit_wave_scheme::rhs_source(coefficients_, grid_cfg_, prev_solution_0, prev_solution_1, source_, boundary_pair,
                                     time, rhs_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution, next_time);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
