#include "lss_ode_implicit_solver_method.hpp"

namespace lss_ode_solvers
{

using d_1d = discretization_1d;

ode_implicit_solver_method::ode_implicit_solver_method(tridiagonal_solver_ptr const &solver_ptr,
                                                       ode_implicit_coefficients_ptr const &coefficients,
                                                       grid_config_1d_ptr const &grid_config)
    : solveru_ptr_{solver_ptr}, coefficients_{coefficients}, grid_cfg_{grid_config}
{
    initialize();
}

void ode_implicit_solver_method::initialize()
{
    low_.resize(coefficients_->space_size_);
    diag_.resize(coefficients_->space_size_);
    high_.resize(coefficients_->space_size_);
    rhs_.resize(coefficients_->space_size_);
}

void ode_implicit_solver_method::split(container_t &low, container_t &diag, container_t &high)
{
    double x{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        x = grid_1d::value(grid_cfg_, t);
        low[t] = coefficients_->A_(x);
        diag[t] = coefficients_->C_(x);
        high[t] = coefficients_->B_(x);
    }
}

ode_implicit_solver_method ::~ode_implicit_solver_method()
{
}

void ode_implicit_solver_method::solve(boundary_1d_pair const &boundary_pair, container_t &solution)
{
    // get the right-hand side of the scheme:
    split(low_, diag_, high_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution);
}

void ode_implicit_solver_method::solve(boundary_1d_pair const &boundary_pair,
                                       std::function<double(double)> const &source, container_t &solution)
{
    // get the right-hand side of the scheme:
    split(low_, diag_, high_);
    d_1d::of_function(grid_cfg_, source, rhs_);
    solveru_ptr_->set_diagonals(low_, diag_, high_);
    solveru_ptr_->set_rhs(rhs_);
    solveru_ptr_->solve(boundary_pair, solution);
}

} // namespace lss_ode_solvers
