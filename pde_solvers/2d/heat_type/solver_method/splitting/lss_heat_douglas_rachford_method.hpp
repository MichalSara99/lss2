/**

    @file      lss_heat_douglas_rachford_method.hpp
    @brief     Douglas-Rachford splitting
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright ? Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_DOUGLAS_RACHFORD_METHOD_HPP_)
#define _LSS_HEAT_DOUGLAS_RACHFORD_METHOD_HPP_

#include <functional>
#include <map>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_enumerations.hpp"
#include "../../../../../containers/lss_matrix_2d.hpp"
#include "../../../../../discretization/lss_grid_config.hpp"
#include "../../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../../implicit_coefficients/lss_heat_coefficients_2d.hpp"
#include "lss_heat_splitting_method.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{
using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::matrix_2d;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;

class douglas_rachford_rhs
{

  public:
    static void rhs_intermed_1(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                               std::size_t const &y_index, double const &y, matrix_2d const &input, double const &time,
                               container_t &solution);

    static void rhs_intermed_1_source(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg,
                                      std::size_t const &y_index, double const &y, matrix_2d const &input,
                                      matrix_2d const &inhom_input, matrix_2d const &inhom_input_next,
                                      double const &time, container_t &solution);

    static void rhs(heat_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_cfg, std::size_t const &x_index,
                    double const &x, matrix_2d const &input, matrix_2d const &inhom_input, double const &time,
                    container_t &solution);
};

/**
    heat_douglas_rachford_method object
 */
class heat_douglas_rachford_method : public heat_splitting_method
{

  private:
    // constants:
    const double cone_ = 1.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solvery_ptr_;
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    heat_coefficients_2d_ptr coefficients_;
    grid_config_2d_ptr grid_cfg_;
    // containers:
    container_t low_, diag_, high_, rhs_;

    explicit heat_douglas_rachford_method() = delete;

    void initialize(bool is_heat_source_set);

    void split_0(double const &y, double const &time, container_t &low, container_t &diag, container_t &high);

    void split_1(double const &x, double const &time, container_t &low, container_t &diag, container_t &high);

  public:
    explicit heat_douglas_rachford_method(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solvery_ptr,
                                          lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru_ptr,
                                          heat_coefficients_2d_ptr const &coefficients,
                                          grid_config_2d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_douglas_rachford_method();

    heat_douglas_rachford_method(heat_douglas_rachford_method const &) = delete;
    heat_douglas_rachford_method(heat_douglas_rachford_method &&) = delete;
    heat_douglas_rachford_method &operator=(heat_douglas_rachford_method const &) = delete;
    heat_douglas_rachford_method &operator=(heat_douglas_rachford_method &&) = delete;

    void solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizontal_boundary_pair,
               boundary_2d_pair const &vertical_boundary_pair, double const &time, matrix_2d &solution) override;

    void solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizontal_boundary_pair,
               boundary_2d_pair const &vertical_boundary_pair, double const &time,
               std::function<double(double, double)> const &heat_source, matrix_2d &solution) override;
};

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_DOUGLAS_RACHFORD_METHOD_HPP_
