/**

    @file      lss_hhw_implicit_boundary_solver.hpp
    @brief     Heston-Hull-White's implicit boundary solver for volatility = 0
    @details   ~
    @author    Michal Sara
    @date      5.03.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HHW_IMPLICIT_BOUNDARY_SOLVER_HPP_)
#define _LSS_HHW_IMPLICIT_BOUNDARY_SOLVER_HPP_

#include <vector>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_enumerations.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../containers/lss_matrix_3d.hpp"
#include "../../../../../discretization/lss_grid_config.hpp"
#include "../../../../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"
#include "../implicit_coefficients/lss_heat_coefficients_3d.hpp"

namespace lss_pde_solvers
{
namespace three_dimensional
{

using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_containers::matrix_3d;
using lss_grids::grid_config_3d_ptr;
using lss_utility::container_t;

/**
    implicit_hhw_boundary_rhs object
 */
class implicit_hhw_boundary_rhs
{

  public:
    static void rhs(heat_coefficients_3d_ptr const &cfg, grid_config_3d_ptr const &grid_cfg, std::size_t const &y_index,
                    double const &y, boundary_3d_pair const &x_boundary_pair, boundary_3d_pair const &z_boundary_pair,
                    matrix_3d const &input, double const &time, matrix_2d &solution);

    static void rhs_source(heat_coefficients_3d_ptr const &cfg, grid_config_3d_ptr const &grid_cfg,
                           std::size_t const &y_index, double const &y, boundary_3d_pair const &x_boundary_pair,
                           boundary_3d_pair const &z_boundary_pair, matrix_3d const &input,
                           matrix_3d const &inhom_input, double const &time, matrix_2d &solution);
};

/**
    hhw_implicit_boundary_solver object
 */
class hhw_implicit_boundary_solver
{

  private:
    // constants:
    const double cone_ = 1.0;
    const double ctwo_ = 2.0;
    // solvers:
    lss_tridiagonal_solver::tridiagonal_solver_ptr solveru_ptr_;
    // scheme coefficients:
    heat_coefficients_3d_ptr coefficients_;
    grid_config_3d_ptr grid_cfg_;
    // containers:
    container_t low_, diag_, high_, rhs_;

    explicit hhw_implicit_boundary_solver() = delete;

  public:
    explicit hhw_implicit_boundary_solver(lss_tridiagonal_solver::tridiagonal_solver_ptr const &solveru,
                                          heat_coefficients_3d_ptr const &coefficients,
                                          grid_config_3d_ptr const &grid_config);

    ~hhw_implicit_boundary_solver();

    hhw_implicit_boundary_solver(hhw_implicit_boundary_solver const &) = delete;
    hhw_implicit_boundary_solver(hhw_implicit_boundary_solver &&) = delete;
    hhw_implicit_boundary_solver &operator=(hhw_implicit_boundary_solver const &) = delete;
    hhw_implicit_boundary_solver &operator=(hhw_implicit_boundary_solver &&) = delete;

    void split(double const &x, double const &y, double const &time, container_t &low, container_t &diag,
               container_t &high);

    void solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
               boundary_3d_ptr const &y_upper_boundary_ptr, boundary_3d_pair const &z_boundary_pair, double const &time,
               matrix_3d &solution);

    void solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
               boundary_3d_pair const &z_boundary_pair, double const &time, matrix_3d &solution);
};

using hhw_implicit_boundary_solver_ptr = sptr_t<hhw_implicit_boundary_solver>;

} // namespace three_dimensional

} // namespace lss_pde_solvers

#endif /// _LSS_HHW_IMPLICIT_BOUNDARY_SOLVER_HPP_
