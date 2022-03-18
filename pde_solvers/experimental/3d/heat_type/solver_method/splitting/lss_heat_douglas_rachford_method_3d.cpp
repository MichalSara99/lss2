#include "lss_heat_douglas_rachford_method_3d.hpp"

#include "../../../../../../common/lss_macros.hpp"
#include "../../../../../../discretization/lss_grid.hpp"
//#include <future>
#include <omp.h>

namespace lss_pde_solvers
{

using lss_grids::grid_3d;

namespace three_dimensional
{

void douglas_rachford_3d_rhs::rhs_intermed_1(heat_coefficients_3d_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                                             std::size_t const &y_index, double const &y, std::size_t const &z_index,
                                             double const &z, matrix_3d const &input, double const &time,
                                             container_t &solution)
{
    auto const one = 1.0;
    auto const two = 2.0;
    auto const &M_1 = cfs->M_1_;
    auto const &M_2 = cfs->M_2_;
    auto const &M_3 = cfs->M_3_;
    auto const &P_1 = cfs->P_1_;
    auto const &P_2 = cfs->P_2_;
    auto const &P_3 = cfs->P_3_;
    auto const &S_1 = cfs->S_1_;
    auto const &S_2 = cfs->S_2_;
    auto const &S_3 = cfs->S_3_;
    auto const &D = cfs->D_;
    auto const &E = cfs->E_;
    auto const &F = cfs->F_;

    auto const beta_1 = cfs->beta_1_;
    auto const beta_2 = cfs->beta_2_;
    auto const beta_3 = cfs->beta_3_;
    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double x{};
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_3d::value_1(grid_cfg, t);
        solution[t] =
            (beta_1 * D(time, x, y, z) *
             (input(t - 1, y_index - 1, z_index) - input(t - 1, y_index + 1, z_index) -
              input(t + 1, y_index - 1, z_index) + input(t + 1, y_index + 1, z_index))) +
            (beta_2 * E(time, x, y, z) *
             (input(t - 1, y_index, z_index - 1) - input(t - 1, y_index, z_index + 1) -
              input(t + 1, y_index, z_index - 1) + input(t + 1, y_index, z_index + 1))) +
            (beta_3 * F(time, x, y, z) *
             (input(t, y_index - 1, z_index - 1) - input(t, y_index - 1, z_index + 1) -
              input(t, y_index + 1, z_index - 1) + input(t, y_index + 1, z_index + 1))) +
            ((one - theta) * M_1(time, x, y, z) * input(t - 1, y_index, z_index)) +
            (M_2(time, x, y, z) * input(t, y_index - 1, z_index)) +
            (M_3(time, x, y, z) * input(t, y_index, z_index - 1)) +
            ((one - theta) * P_1(time, x, y, z) * input(t + 1, y_index, z_index)) +
            (P_2(time, x, y, z) * input(t, y_index + 1, z_index)) +
            (P_3(time, x, y, z) * input(t, y_index, z_index + 1)) +
            ((one - two * S_2(time, x, y, z) - two * S_3(time, x, y, z) - two * (one - theta) * S_1(time, x, y, z)) *
             input(t, y_index, z_index));
    }
}

void douglas_rachford_3d_rhs::rhs_intermed_1_source(heat_coefficients_3d_ptr const &cfs,
                                                    grid_config_3d_ptr const &grid_cfg, std::size_t const &y_index,
                                                    double const &y, std::size_t const &z_index, double const &z,
                                                    matrix_3d const &input, matrix_3d const &inhom_input,
                                                    matrix_3d const &inhom_input_next, double const &time,
                                                    container_t &solution)
{
    auto const one = 1.0;
    auto const two = 2.0;
    auto const &M_1 = cfs->M_1_;
    auto const &M_2 = cfs->M_2_;
    auto const &M_3 = cfs->M_3_;
    auto const &P_1 = cfs->P_1_;
    auto const &P_2 = cfs->P_2_;
    auto const &P_3 = cfs->P_3_;
    auto const &S_1 = cfs->S_1_;
    auto const &S_2 = cfs->S_2_;
    auto const &S_3 = cfs->S_3_;
    auto const &D = cfs->D_;
    auto const &E = cfs->E_;
    auto const &F = cfs->F_;

    auto const beta_1 = cfs->beta_1_;
    auto const beta_2 = cfs->beta_2_;
    auto const beta_3 = cfs->beta_3_;
    auto const rho = cfs->rho_;
    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double x{};
    for (std::size_t t = 1; t < N; ++t)
    {
        x = grid_3d::value_1(grid_cfg, t);
        solution[t] =
            (beta_1 * D(time, x, y, z) *
             (input(t - 1, y_index - 1, z_index) - input(t - 1, y_index + 1, z_index) -
              input(t + 1, y_index - 1, z_index) + input(t + 1, y_index + 1, z_index))) +
            (beta_2 * E(time, x, y, z) *
             (input(t - 1, y_index, z_index - 1) - input(t - 1, y_index, z_index + 1) -
              input(t + 1, y_index, z_index - 1) + input(t + 1, y_index, z_index + 1))) +
            (beta_3 * F(time, x, y, z) *
             (input(t, y_index - 1, z_index - 1) - input(t, y_index - 1, z_index + 1) -
              input(t, y_index + 1, z_index - 1) + input(t, y_index + 1, z_index + 1))) +
            ((one - theta) * M_1(time, x, y, z) * input(t - 1, y_index, z_index)) +
            (M_2(time, x, y, z) * input(t, y_index - 1, z_index)) +
            (M_3(time, x, y, z) * input(t, y_index, z_index - 1)) +
            ((one - theta) * P_1(time, x, y, z) * input(t + 1, y_index, z_index)) +
            (P_2(time, x, y, z) * input(t, y_index + 1, z_index)) +
            (P_3(time, x, y, z) * input(t, y_index, z_index + 1)) +
            ((one - two * S_2(time, x, y, z) - two * S_3(time, x, y, z) - two * (one - theta) * S_1(time, x, y, z)) *
             input(t, y_index, z_index)) +
            rho * (one - theta) * inhom_input(t, y_index, z_index) +
            rho * theta * inhom_input_next(t, y_index, z_index);
    }
}

void douglas_rachford_3d_rhs::rhs_intermed_2(heat_coefficients_3d_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                                             std::size_t const &x_index, double const &x, std::size_t const &z_index,
                                             double const &z, matrix_3d const &input, matrix_3d const &inhom_input,
                                             double const &time, container_t &solution)
{
    auto const two = 2.0;
    auto const &M_2 = cfs->M_2_;
    auto const &P_2 = cfs->P_2_;
    auto const &S_2 = cfs->S_2_;
    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double y{};
    for (std::size_t t = 1; t < N; ++t)
    {
        y = grid_3d::value_2(grid_cfg, t);
        solution[t] = (-theta * M_2(time, x, y, z) * input(x_index, t - 1, z_index)) +
                      (theta * two * S_2(time, x, y, z) * input(x_index, t, z_index)) -
                      (theta * P_2(time, x, y, z) * input(x_index, t + 1, z_index)) + inhom_input(x_index, t, z_index);
    }
}

void douglas_rachford_3d_rhs::rhs(heat_coefficients_3d_ptr const &cfs, grid_config_3d_ptr const &grid_cfg,
                                  std::size_t const &x_index, double const &x, std::size_t const &y_index,
                                  double const &y, matrix_3d const &input, matrix_3d const &inhom_input,
                                  double const &time, container_t &solution)
{
    auto const two = 2.0;
    auto const &M_3 = cfs->M_2_;
    auto const &P_3 = cfs->P_2_;
    auto const &S_3 = cfs->S_2_;
    auto const theta = cfs->theta_;

    const std::size_t N = solution.size() - 1;
    double z{};
    for (std::size_t t = 1; t < N; ++t)
    {
        z = grid_3d::value_3(grid_cfg, t);
        solution[t] = (-theta * M_3(time, x, y, z) * input(x_index, y_index, t - 1)) +
                      (theta * two * S_3(time, x, y, z) * input(x_index, y_index, t)) -
                      (theta * P_3(time, x, y, z) * input(x_index, y_index, t + 1)) + inhom_input(x_index, y_index, t);
    }
}

void heat_douglas_rachford_method_3d::split_0(double const &y, double const &z, double const &time, container_t &low,
                                              container_t &diag, container_t &high)
{
    double x{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        x = grid_3d::value_1(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->M_1_(time, x, y, z));
        diag[t] = (cone_ + ctwo_ * coefficients_->theta_ * coefficients_->S_1_(time, x, y, z));
        high[t] = (-coefficients_->theta_ * coefficients_->P_1_(time, x, y, z));
    }
}

void heat_douglas_rachford_method_3d::split_1(double const &x, double const &z, double const &time, container_t &low,
                                              container_t &diag, container_t &high)
{
    double y{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        y = grid_3d::value_2(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->M_2_(time, x, y, z));
        diag[t] = (cone_ + ctwo_ * coefficients_->theta_ * coefficients_->S_2_(time, x, y, z));
        high[t] = (-coefficients_->theta_ * coefficients_->P_2_(time, x, y, z));
    }
}

void heat_douglas_rachford_method_3d::split_2(double const &x, double const &y, double const &time, container_t &low,
                                              container_t &diag, container_t &high)
{
    double z{};
    for (std::size_t t = 0; t < low.size(); ++t)
    {
        z = grid_3d::value_3(grid_cfg_, t);
        low[t] = (-coefficients_->theta_ * coefficients_->M_3_(time, x, y, z));
        diag[t] = (cone_ + ctwo_ * coefficients_->theta_ * coefficients_->S_3_(time, x, y, z));
        high[t] = (-coefficients_->theta_ * coefficients_->P_3_(time, x, y, z));
    }
}

heat_douglas_rachford_method_3d::heat_douglas_rachford_method_3d(
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery1_ptr,
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery2_ptr,
    lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru_ptr, heat_coefficients_3d_ptr const &coefficients,
    grid_config_3d_ptr const &grid_config, bool is_heat_source_set)
    : solvery1_ptr_{solvery1_ptr}, solvery2_ptr_{solvery2_ptr}, solveru_ptr_{solveru_ptr},
      coefficients_{coefficients}, grid_cfg_{grid_config}
{
    initialize(is_heat_source_set);
}

heat_douglas_rachford_method_3d::~heat_douglas_rachford_method_3d()
{
}

void heat_douglas_rachford_method_3d::initialize(bool is_heat_source_set)
{
}

void heat_douglas_rachford_method_3d::solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
                                            boundary_3d_pair const &y_boundary_pair,
                                            boundary_3d_pair const &z_boundary_pair, double const &time,
                                            matrix_3d &solution)
{
    // 3D container for intermediate solution:
    matrix_3d inter_solution_1(coefficients_->space_size_x_, coefficients_->space_size_y_, coefficients_->space_size_z_,
                               double{});
    // 1D container for intermediate solution:
    container_t solution_v(double{}, coefficients_->space_size_x_);
    low_.resize(coefficients_->space_size_x_);
    diag_.resize(coefficients_->space_size_x_);
    high_.resize(coefficients_->space_size_x_);
    rhs_.resize(coefficients_->space_size_x_);
    double x{}, y{}, z{};
    const std::size_t total_yz = coefficients_->space_size_y_ * coefficients_->space_size_z_;
    //  TODO: this may be paralelised ?
    for (std::size_t t = 0, j = 0, k = 0; t < total_yz; ++t)
    {
        j = t / coefficients_->space_size_z_;
        k = t - j * coefficients_->space_size_z_;
        if (j == 0 || k == 0)
            continue;
        if (j == (coefficients_->space_size_y_ - 1) || k == (coefficients_->space_size_z_ - 1))
            continue;
        y = grid_3d::value_2(grid_cfg_, j);
        z = grid_3d::value_3(grid_cfg_, k);
        split_0(y, z, time, low_, diag_, high_);
        douglas_rachford_3d_rhs::rhs_intermed_1(coefficients_, grid_cfg_, j, y, k, z, prev_solution, time, rhs_);
        solvery1_ptr_->set_diagonals(low_, diag_, high_);
        solvery1_ptr_->set_rhs(rhs_);
        solvery1_ptr_->solve(x_boundary_pair, solution_v, time, y, z);
        inter_solution_1.column_array(j, k) = solution_v;
    }
    // 3D container for intermediate solution:
    matrix_3d inter_solution_2(coefficients_->space_size_x_, coefficients_->space_size_y_, coefficients_->space_size_z_,
                               double{});
    // 1D container for intermediate solution:
    solution_v.resize(coefficients_->space_size_y_);
    // containers for second split solver:
    low_.resize(coefficients_->space_size_y_);
    diag_.resize(coefficients_->space_size_y_);
    high_.resize(coefficients_->space_size_y_);
    rhs_.resize(coefficients_->space_size_y_);
    const std::size_t total_xz = coefficients_->space_size_x_ * coefficients_->space_size_z_;
    for (std::size_t t = 0, i = 0, k = 0; t < total_xz; ++t)
    {
        i = t / coefficients_->space_size_z_;
        k = t - i * coefficients_->space_size_z_;
        if (i == 0 || k == 0)
            continue;
        if (i == (coefficients_->space_size_x_ - 1) || k == (coefficients_->space_size_z_ - 1))
            continue;
        x = grid_3d::value_1(grid_cfg_, i);
        z = grid_3d::value_3(grid_cfg_, k);
        split_1(x, z, time, low_, diag_, high_);
        douglas_rachford_3d_rhs::rhs_intermed_2(coefficients_, grid_cfg_, i, x, k, z, prev_solution, inter_solution_1,
                                                time, rhs_);
        solvery2_ptr_->set_diagonals(low_, diag_, high_);
        solvery2_ptr_->set_rhs(rhs_);
        solvery2_ptr_->solve(y_boundary_pair, solution_v, time, x, z);
        inter_solution_2.row_array(i, k) = solution_v;
    }

    // 1D container for final solution:
    solution_v.resize(coefficients_->space_size_z_);
    // containers for second split solver:
    low_.resize(coefficients_->space_size_z_);
    diag_.resize(coefficients_->space_size_z_);
    high_.resize(coefficients_->space_size_z_);
    rhs_.resize(coefficients_->space_size_z_);
    const std::size_t total_xy = coefficients_->space_size_x_ * coefficients_->space_size_y_;
    // TODO: this may be paralelised ?
    for (std::size_t t = 0, i = 0, j = 0; t < total_xy; ++t)
    {
        i = t / coefficients_->space_size_y_;
        j = t - i * coefficients_->space_size_y_;
        if (i == 0 || j == 0)
            continue;
        if (i == (coefficients_->space_size_x_ - 1) || j == (coefficients_->space_size_y_ - 1))
            continue;
        x = grid_3d::value_1(grid_cfg_, i);
        y = grid_3d::value_2(grid_cfg_, j);
        split_2(x, y, time, low_, diag_, high_);
        douglas_rachford_3d_rhs::rhs(coefficients_, grid_cfg_, i, x, j, y, prev_solution, inter_solution_2, time, rhs_);
        solveru_ptr_->set_diagonals(low_, diag_, high_);
        solveru_ptr_->set_rhs(rhs_);
        solveru_ptr_->solve(z_boundary_pair, solution_v, time, x, y);
        solution.layer_array(i, j) = solution_v;
    }
}

void heat_douglas_rachford_method_3d::solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
                                            boundary_3d_pair const &y_boundary_pair,
                                            boundary_3d_pair const &z_boundary_pair, double const &time,
                                            std::function<double(double, double, double)> const &heat_source,
                                            matrix_3d &solution)
{
}

} // namespace three_dimensional

} // namespace lss_pde_solvers
