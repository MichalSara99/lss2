#include "lss_tridiagonal_solver.hpp"

namespace lss_tridiagonal_solver
{

tridiagonal_solver::tridiagonal_solver(std::size_t discretization_size, factorization_enum factorization)
    : discretization_size_{discretization_size}, factorization_{factorization}
{
}

void tridiagonal_solver::set_diagonals(container_t lower_diagonal, container_t diagonal, container_t upper_diagonal)
{
    LSS_ASSERT(lower_diagonal.size() == discretization_size_, "Inncorect size for lower_diagonal");
    LSS_ASSERT(diagonal.size() == discretization_size_, "Inncorect size for diagonal");
    LSS_ASSERT(upper_diagonal.size() == discretization_size_, "Inncorect size for upper_diagonal");
    a_ = std::move(lower_diagonal);
    b_ = std::move(diagonal);
    c_ = std::move(upper_diagonal);
}

void tridiagonal_solver::set_rhs(container_t const &rhs)
{
    LSS_ASSERT(rhs.size() == discretization_size_, "Inncorect size for right-hand side");
    f_ = rhs;
}

void tridiagonal_solver::set_factorization(factorization_enum factorization)
{
    factorization_ = factorization;
}

void tridiagonal_solver::set_omega(double value)
{
    omega_ = value;
}

void tridiagonal_solver::solve(boundary_1d_pair const &boundary, container_t &solution)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, solution, factorization_, double{});
}

void tridiagonal_solver::solve(boundary_1d_pair const &boundary, container_t &solution, double at_time)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, solution, factorization_, at_time);
}

void tridiagonal_solver::solve(boundary_2d_pair const &boundary, container_t &solution, double at_time,
                               double space_arg)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, solution, factorization_, at_time, space_arg);
}

void tridiagonal_solver::solve(boundary_3d_pair const &boundary, container_t &solution, double at_time,
                               double space_1_arg, double space_2_arg)
{
    LSS_ASSERT(solution.size() == discretization_size_, "Incorrect size of solution container");
    kernel(boundary, solution, factorization_, at_time, space_1_arg, space_2_arg);
}

} // namespace lss_tridiagonal_solver
