#include "lss_heat_euler_solver_method_2d.hpp"

#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

using lss_grids::grid_2d;
using d_2d = discretization_2d;

namespace two_dimensional
{

void explicit_euler_rhs::rhs(heat_euler_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                             matrix_2d const &input, double const &time, matrix_2d &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P = cfs->P_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &Z = cfs->Z_;
    auto const &W = cfs->W_;
    auto const &C = cfs->C_;
    auto const rows = input.rows() - 1;
    auto const cols = input.columns() - 1;

    double x{}, y{}, val{};
    for (std::size_t r = 1; r < rows; ++r)
    {
        x = grid_2d::value_1(grid_config, r);
        for (std::size_t c = 1; c < cols; ++c)
        {
            y = grid_2d::value_2(grid_config, c);
            val = C(time, x, y) * input(r - 1, c - 1) + M(time, x, y) * input(r - 1, c) -
                  C(time, x, y) * input(r - 1, c + 1) + M_tilde(time, x, y) * input(r, c - 1) +
                  (one - Z(time, x, y) - W(time, x, y)) * input(r, c) + P_tilde(time, x, y) * input(r, c + 1) -
                  C(time, x, y) * input(r + 1, c - 1) + P(time, x, y) * input(r + 1, c) +
                  C(time, x, y) * input(r + 1, c + 1);
            solution(r, c, val);
        }
    }
}

void explicit_euler_rhs::rhs_source(heat_euler_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                                    matrix_2d const &input, double const &time, matrix_2d const &inhom_input,
                                    matrix_2d &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P = cfs->P_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &Z = cfs->Z_;
    auto const &W = cfs->W_;
    auto const &C = cfs->C_;
    auto const rho = cfs->rho_;
    auto const rows = input.rows() - 1; // without boundary
    auto const cols = input.columns() - 1;

    double x{}, y{}, val{};
    for (std::size_t r = 1; r < rows; ++r)
    {
        x = grid_2d::value_1(grid_config, r);
        for (std::size_t c = 1; c < cols; ++c)
        {
            y = grid_2d::value_2(grid_config, c);
            val = C(time, x, y) * input(r - 1, c - 1) + M(time, x, y) * input(r - 1, c) -
                  C(time, x, y) * input(r - 1, c + 1) + M_tilde(time, x, y) * input(r, c - 1) +
                  (one - Z(time, x, y) - W(time, x, y)) * input(r, c) + P_tilde(time, x, y) * input(r, c + 1) -
                  C(time, x, y) * input(r + 1, c - 1) + P(time, x, y) * input(r + 1, c) +
                  C(time, x, y) * input(r + 1, c + 1) + rho * inhom_input(r, c);
            solution(r, c, val);
        }
    }
}

void heat_euler_solver_method_2d::initialize(bool is_heat_source_set)
{
    if (is_heat_source_set)
    {
        source_ = std::make_shared<matrix_2d>(coefficients_->space_size_x_, coefficients_->space_size_y_);
    }
}

heat_euler_solver_method_2d::heat_euler_solver_method_2d(heat_euler_coefficients_2d_ptr const &coefficients,
                                                         grid_config_2d_ptr const &grid_config, bool is_heat_source_set)
    : coefficients_{coefficients}, heat_explicit_solver_method_2d(grid_config)
{
    initialize(is_heat_source_set);
}

heat_euler_solver_method_2d::~heat_euler_solver_method_2d()
{
}

void heat_euler_solver_method_2d::solve(matrix_2d &prev_solution, double const &time, matrix_2d &solution)
{
    explicit_euler_rhs::rhs(coefficients_, grid_cfg_, prev_solution, time, solution);
}

void heat_euler_solver_method_2d::solve(matrix_2d &prev_solution, double const &time,
                                        std::function<double(double, double, double)> const &heat_source,
                                        matrix_2d &solution)
{
    d_2d::of_function(grid_cfg_, time, heat_source, *source_);
    explicit_euler_rhs::rhs_source(coefficients_, grid_cfg_, prev_solution, time, *source_, solution);
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
