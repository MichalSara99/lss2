#include "lss_heat_douglas_rachford_method.hpp"

#include "../../../../common/lss_macros.hpp"
#include "../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

using lss_grids::grid_2d;

namespace two_dimensional
{

void implicit_heston_scheme::rhs_intermed_1(heston_implicit_coefficients_ptr const &cfs,
                                            grid_config_2d_ptr const &grid_cfg, std::size_t const &y_index,
                                            double const &y, container_2d<by_enum::Row> const &input,
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

void implicit_heston_scheme::rhs_intermed_1_source(heston_implicit_coefficients_ptr const &cfs,
                                                   grid_config_2d_ptr const &grid_cfg, std::size_t const &y_index,
                                                   double const &y, container_2d<by_enum::Row> const &input,
                                                   container_2d<by_enum::Row> const &inhom_input,
                                                   container_2d<by_enum::Row> const &inhom_input_next,
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

void implicit_heston_scheme::rhs(heston_implicit_coefficients_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                 std::size_t const &x_index, double const &x, container_2d<by_enum::Row> const &input,
                                 container_2d<by_enum::Row> const &inhom_input, double const &time,
                                 container_t &solution)
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

void heat_douglas_rachford_method::initialize(bool is_heat_source_set)
{
}

void heat_douglas_rachford_method::split_0(double const &y, double const &time, container_t &low, container_t &diag,
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

void heat_douglas_rachford_method::split_1(double const &x, double const &time, container_t &low, container_t &diag,
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

heat_douglas_rachford_method::heat_douglas_rachford_method(
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery_ptr,
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru_ptr,
    heston_implicit_coefficients_ptr const &coefficients, grid_config_2d_ptr const &grid_config,
    bool is_heat_source_set)
    : solvery_ptr_{solvery_ptr}, solveru_ptr_{solveru_ptr}, coefficients_{coefficients}, grid_cfg_{grid_config}
{
    initialize(is_heat_source_set);
}

heat_douglas_rachford_method::~heat_douglas_rachford_method()
{
}

void heat_douglas_rachford_method::solve(container_2d<by_enum::Row> const &prev_solution,
                                         boundary_2d_pair const &horizontal_boundary_pair,
                                         boundary_2d_pair const &vertical_boundary_pair, double const &time,
                                         container_2d<by_enum::Row> &solution)
{
    // 2D container for intermediate solution:
    container_2d<by_enum::Column> inter_solution(coefficients_->space_size_x_, coefficients_->space_size_y_, double{});
    // 1D container for intermediate solution:
    container_t solution_v(coefficients_->space_size_x_, double{});
    low_.resize(coefficients_->space_size_x_);
    diag_.resize(coefficients_->space_size_x_);
    high_.resize(coefficients_->space_size_x_);
    rhs_.resize(coefficients_->space_size_x_);
    double y{};
    for (std::size_t j = 1; j < coefficients_->space_size_y_ - 1; ++j)
    {
        y = grid_2d::value_2(grid_cfg_, j);
        split_0(y, time, low_, diag_, high_);
        implicit_heston_scheme::rhs_intermed_1(coefficients_, grid_cfg_, j, y, prev_solution, time, rhs_);
        solvery_ptr_->set_diagonals(low_, diag_, high_);
        solvery_ptr_->set_rhs(rhs_);
        solvery_ptr_->solve(horizontal_boundary_pair, solution_v, time, y);
        inter_solution(j, solution_v);
    }

    // 1D container for final solution:
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
        implicit_heston_scheme::rhs(coefficients_, grid_cfg_, i, x, prev_solution, inter_solution, time, rhs_);
        solveru_ptr_->set_diagonals(low_, diag_, high_);
        solveru_ptr_->set_rhs(rhs_);
        solveru_ptr_->solve(vertical_boundary_pair, solution_v, time, x);
        solution(i, solution_v);
    }
}

void heat_douglas_rachford_method::solve(container_2d<by_enum::Row> const &prev_solution,
                                         boundary_2d_pair const &horizontal_boundary_pair,
                                         boundary_2d_pair const &vertical_boundary_pair, double const &time,
                                         std::function<double(double, double)> const &heat_source,
                                         container_2d<by_enum::Row> &solution)
{
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
