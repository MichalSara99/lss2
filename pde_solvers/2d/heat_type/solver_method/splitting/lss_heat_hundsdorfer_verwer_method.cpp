#include "lss_heat_hundsdorfer_verwer_method.hpp"

#include "../../../../../common/lss_macros.hpp"
#include "../../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

using lss_grids::grid_2d;

namespace two_dimensional
{

void hundsdorfer_verwer_rhs::rhs_intermed_1(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                            std::size_t const &y_index, double const &y, matrix_2d const &input,
                                            double const &time, container_t &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P = cfs->P_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &Z = cfs->Z_;
    auto const &W = cfs->W_;
    auto const &C = cfs->C_;
    auto gamma = cfs->gamma_;

    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double x{};
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_2d::value_1(grid_cfg, t);
        solution[t] =
            (gamma * C(time, x, y) * input(t - 1, y_index - 1)) +
            ((one - theta) * M(time, x, y) * input(t - 1, y_index)) -
            (gamma * C(time, x, y) * input(t - 1, y_index + 1)) + (M_tilde(time, x, y) * input(t, y_index - 1)) +
            ((one - W(time, x, y) - (one - theta) * Z(time, x, y)) * input(t, y_index)) +
            (P_tilde(time, x, y) * input(t, y_index + 1)) - (gamma * C(time, x, y) * input(t + 1, y_index - 1)) +
            ((one - theta) * P(time, x, y) * input(t + 1, y_index)) +
            (gamma * C(time, x, y) * input(t + 1, y_index + 1));
    }
}

void hundsdorfer_verwer_rhs::rhs_intermed_1_source(heat_coefficients_2d_ptr const &cfs,
                                                   grid_config_2d_ptr const &grid_cfg, std::size_t const &y_index,
                                                   double const &y, matrix_2d const &input,
                                                   matrix_2d const &inhom_input, matrix_2d const &inhom_input_next,
                                                   double const &time, container_t &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P = cfs->P_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &Z = cfs->Z_;
    auto const &W = cfs->W_;
    auto const &C = cfs->C_;

    auto const gamma = cfs->gamma_;
    auto const theta = cfs->theta_;
    auto const rho = cfs->rho_;

    const std::size_t N = solution.size() - 1;
    double x{};
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_2d::value_1(grid_cfg, t);
        solution[t] =
            (gamma * C(time, x, y) * input(t - 1, y_index - 1)) +
            ((one - theta) * M(time, x, y) * input(t - 1, y_index)) -
            (gamma * C(time, x, y) * input(t - 1, y_index + 1)) + (M_tilde(time, x, y) * input(t, y_index - 1)) +
            ((one - W(time, x, y) - (one - theta) * Z(time, x, y)) * input(t, y_index)) +
            (P_tilde(time, x, y) * input(t, y_index + 1)) - (gamma * C(time, x, y) * input(t + 1, y_index - 1)) +
            ((one - theta) * P(time, x, y) * input(t + 1, y_index)) +
            (gamma * C(time, x, y) * input(t + 1, y_index + 1)) + (theta * rho * inhom_input_next(t, y_index)) +
            ((one - theta) * rho * inhom_input(t, y_index));
    }
}

void hundsdorfer_verwer_rhs::rhs_intermed_2(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                            std::size_t const &x_index, double const &x, matrix_2d const &input,
                                            matrix_2d const &inhom_input, double const &time, container_t &solution)
{
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &W = cfs->W_;

    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double y{};
    for (std::size_t t = 1; t < N; ++t)
    {
        y = grid_2d::value_2(grid_cfg, t);
        solution[t] = (-theta * M_tilde(time, x, y) * input(x_index, t - 1)) +
                      (theta * W(time, x, y) * input(x_index, t)) -
                      (theta * P_tilde(time, x, y) * input(x_index, t + 1)) + inhom_input(x_index, t);
    }
}

void hundsdorfer_verwer_rhs::rhs_intermed_3(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                            matrix_2d const &input, matrix_2d const &inhom_input,
                                            matrix_2d const &inhom_input_next, double const &time, matrix_2d &solution)
{
    auto const one = 1.0;
    auto const &M = cfs->M_;
    auto const &P = cfs->P_;
    auto const &Z = cfs->Z_;
    auto const &C = cfs->C_;
    auto const &W = cfs->W_;
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P_tilde = cfs->P_tilde_;

    auto const gamma = cfs->gamma_;
    auto const theta = cfs->theta_;
    auto const zeta = cfs->zeta_;

    const std::size_t rows = solution.rows() - 1;
    const std::size_t cols = solution.columns() - 1;
    double x{}, y{}, val{};
    for (std::size_t r = 1; r < rows; ++r)
    {
        x = grid_2d::value_1(grid_cfg, r);
        for (std::size_t c = 1; c < cols; ++c)
        {
            y = grid_2d::value_2(grid_cfg, c);
            val = (-theta * M(time, x, y) * inhom_input(r - 1, c)) +
                  ((one + theta * Z(time, x, y)) * inhom_input(r, c)) -
                  (theta * P(time, x, y) * inhom_input(r + 1, c)) + (theta * M(time, x, y) * input(r - 1, c)) -
                  (theta * Z(time, x, y) * input(r, c)) + (theta * P(time, x, y) * input(r + 1, c)) +
                  (zeta * gamma * C(time, x, y) *
                   ((inhom_input_next(r + 1, c + 1) - input(r + 1, c + 1)) -
                    (inhom_input_next(r + 1, c - 1) - input(r + 1, c - 1)) -
                    (inhom_input_next(r - 1, c + 1) - input(r - 1, c + 1)) +
                    (inhom_input_next(r - 1, c - 1) - input(r - 1, c - 1)))) +
                  (zeta * (M(time, x, y) * (inhom_input_next(r - 1, c) - input(r - 1, c)) -
                           (W(time, x, y) + Z(time, x, y)) * (inhom_input_next(r, c) - input(r, c)) +
                           P(time, x, y) * (inhom_input_next(r + 1, c) - input(r + 1, c)) +
                           M_tilde(time, x, y) * (inhom_input_next(r, c - 1) - input(r, c - 1)) +
                           P_tilde(time, x, y) * (inhom_input_next(r, c + 1) - input(r, c + 1))));
            solution(r, c, val);
        }
    }
}

void hundsdorfer_verwer_rhs::rhs_intermed_4(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                            std::size_t const &y_index, double const &y, matrix_2d const &input,
                                            matrix_2d const &inhom_input, double const &time, container_t &solution)
{
    auto const &M = cfs->M_;
    auto const &P = cfs->P_;
    auto const &Z = cfs->Z_;

    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double x{};
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_2d::value_1(grid_cfg, t);
        solution[t] = (-theta * M(time, x, y) * input(t - 1, y_index)) + (theta * Z(time, x, y) * input(t, y_index)) -
                      (theta * P(time, x, y) * input(t + 1, y_index)) + inhom_input(t, y_index);
    }
}

void hundsdorfer_verwer_rhs::rhs(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                 std::size_t const &x_index, double const &x, matrix_2d const &input,
                                 matrix_2d const &inhom_input, double const &time, container_t &solution)
{
    auto const &M_tilde = cfs->M_tilde_;
    auto const &P_tilde = cfs->P_tilde_;
    auto const &W = cfs->W_;

    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double y{};
    for (std::size_t t = 1; t < N; ++t)
    {
        y = grid_2d::value_2(grid_cfg, t);
        solution[t] = (-theta * M_tilde(time, x, y) * input(x_index, t - 1)) +
                      (theta * W(time, x, y) * input(x_index, t)) -
                      (theta * P_tilde(time, x, y) * input(x_index, t + 1)) + inhom_input(x_index, t);
    }
}

void heat_hundsdorfer_verwer_method::initialize(bool is_heat_source_set)
{
}

void heat_hundsdorfer_verwer_method::split_0(double const &y, double const &time, container_t &low, container_t &diag,
                                             container_t &high)
{
    double x{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        x = grid_2d::value_1(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->M_(time, x, y));
        diag[t] = (cone_ + coefficients_->theta_ * coefficients_->Z_(time, x, y));
        high[t] = (-coefficients_->theta_ * coefficients_->P_(time, x, y));
    }
}

void heat_hundsdorfer_verwer_method::split_1(double const &x, double const &time, container_t &low, container_t &diag,
                                             container_t &high)
{
    double y{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        y = grid_2d::value_2(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->M_tilde_(time, x, y));
        diag[t] = (cone_ + coefficients_->theta_ * coefficients_->W_(time, x, y));
        high[t] = (-coefficients_->theta_ * coefficients_->P_tilde_(time, x, y));
    }
}

heat_hundsdorfer_verwer_method::heat_hundsdorfer_verwer_method(
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery_ptr,
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru_ptr, heat_coefficients_2d_ptr const &coefficients,
    grid_config_2d_ptr const &grid_config, bool is_heat_source_set)
    : solvery_ptr_{solvery_ptr}, solveru_ptr_{solveru_ptr}, coefficients_{coefficients}, grid_cfg_{grid_config}
{
    initialize(is_heat_source_set);
}

heat_hundsdorfer_verwer_method::~heat_hundsdorfer_verwer_method()
{
}

void heat_hundsdorfer_verwer_method::solve(matrix_2d const &prev_solution,
                                           boundary_2d_pair const &horizontal_boundary_pair,
                                           boundary_2d_pair const &vertical_boundary_pair, double const &time,
                                           matrix_2d &solution)
{
    // 2D container for intermediate solution Y_1:
    matrix_2d inter_solution_1(coefficients_->space_size_x_, coefficients_->space_size_y_, double{});
    // 1D container for intermediate solution Y_1:
    container_t solution_v(double{}, coefficients_->space_size_x_);
    // containers for first split solver:
    low_.resize(coefficients_->space_size_x_);
    diag_.resize(coefficients_->space_size_x_);
    high_.resize(coefficients_->space_size_x_);
    rhs_.resize(coefficients_->space_size_x_);
    double y{};
    for (std::size_t j = 1; j < coefficients_->space_size_y_ - 1; ++j)
    {
        y = grid_2d::value_2(grid_cfg_, j);
        split_0(y, time, low_, diag_, high_);
        hundsdorfer_verwer_rhs::rhs_intermed_1(coefficients_, grid_cfg_, j, y, prev_solution, time, rhs_);
        solvery_ptr_->set_diagonals(low_, diag_, high_);
        solvery_ptr_->set_rhs(rhs_);
        solvery_ptr_->solve(horizontal_boundary_pair, solution_v, time, y);
        inter_solution_1.column(j) = solution_v;
    }

    // 2D container for intermediate solution Y_2:
    matrix_2d inter_solution_2(coefficients_->space_size_x_, coefficients_->space_size_y_, double{});
    // 1D container for intermediate solution Y_2:
    solution_v.resize(coefficients_->space_size_y_);
    // containers for second split solver:
    low_.resize(coefficients_->space_size_y_);
    diag_.resize(coefficients_->space_size_y_);
    high_.resize(coefficients_->space_size_y_);
    rhs_.resize(coefficients_->space_size_y_);
    double x{};
    for (std::size_t i = 1; i < coefficients_->space_size_x_ - 1; ++i)
    {
        x = grid_2d::value_1(grid_cfg_, i);
        split_1(x, time, low_, diag_, high_);
        hundsdorfer_verwer_rhs::rhs_intermed_2(coefficients_, grid_cfg_, i, x, prev_solution, inter_solution_1, time,
                                               rhs_);
        solveru_ptr_->set_diagonals(low_, diag_, high_);
        solveru_ptr_->set_rhs(rhs_);
        solveru_ptr_->solve(vertical_boundary_pair, solution_v, time, x);
        inter_solution_2.row(i) = solution_v;
    }

    // 2D container for intermediate solution Y_3:
    matrix_2d inter_solution_3(coefficients_->space_size_x_, coefficients_->space_size_y_, double{});
    hundsdorfer_verwer_rhs::rhs_intermed_3(coefficients_, grid_cfg_, prev_solution, inter_solution_1, inter_solution_2,
                                           time, inter_solution_3);

    // 1D container for intermediate solution Y_4:
    solution_v.resize(coefficients_->space_size_x_);
    // containers for second split solver:
    low_.resize(coefficients_->space_size_x_);
    diag_.resize(coefficients_->space_size_x_);
    high_.resize(coefficients_->space_size_x_);
    rhs_.resize(coefficients_->space_size_x_);
    for (std::size_t j = 1; j < coefficients_->space_size_y_ - 1; ++j)
    {
        y = grid_2d::value_2(grid_cfg_, j);
        split_0(y, time, low_, diag_, high_);
        hundsdorfer_verwer_rhs::rhs_intermed_4(coefficients_, grid_cfg_, j, y, inter_solution_2, inter_solution_3, time,
                                               rhs_);
        solvery_ptr_->set_diagonals(low_, diag_, high_);
        solvery_ptr_->set_rhs(rhs_);
        solvery_ptr_->solve(horizontal_boundary_pair, solution_v, time, y);
        inter_solution_1.column(j) = solution_v;
    }

    // 1D container for final solution:
    solution_v.resize(coefficients_->space_size_y_);
    // containers for second split solver:
    low_.resize(coefficients_->space_size_y_);
    diag_.resize(coefficients_->space_size_y_);
    high_.resize(coefficients_->space_size_y_);
    rhs_.resize(coefficients_->space_size_y_);
    for (std::size_t i = 1; i < coefficients_->space_size_x_ - 1; ++i)
    {
        x = grid_2d::value_1(grid_cfg_, i);
        split_1(x, time, low_, diag_, high_);
        hundsdorfer_verwer_rhs::rhs(coefficients_, grid_cfg_, i, x, inter_solution_2, inter_solution_1, time, rhs_);
        solveru_ptr_->set_diagonals(low_, diag_, high_);
        solveru_ptr_->set_rhs(rhs_);
        solveru_ptr_->solve(vertical_boundary_pair, solution_v, time, x);
        solution.row(i) = solution_v;
    }
}

void heat_hundsdorfer_verwer_method::solve(matrix_2d const &prev_solution,
                                           boundary_2d_pair const &horizontal_boundary_pair,
                                           boundary_2d_pair const &vertical_boundary_pair, double const &time,
                                           std::function<double(double, double)> const &heat_source,
                                           matrix_2d &solution)
{
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
