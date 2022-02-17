/**

    @file      lss_heston_boundary_solver.hpp
    @brief     Heston's boundary solver for volatility = 0
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HESTON_BOUNDARY_SOLVER_HPP_)
#define _LSS_HESTON_BOUNDARY_SOLVER_HPP_

#include <vector>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../implicit_coefficients/lss_heat_coefficients_2d.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::matrix_2d;
using lss_grids::grid_config_2d_ptr;
using lss_utility::container_t;

/**
    heston_boundary_rhs object
 */
class heston_boundary_rhs
{

  public:
    static void rhs(heat_coefficients_2d_ptr const &cfg, grid_config_2d_ptr const &grid_cfg, std::size_t const &y_index,
                    double const &y, boundary_2d_pair const &horizontal_boundary_pair, matrix_2d const &input,
                    double const &time, container_t &solution);

    static void rhs_source(heat_coefficients_2d_ptr const &cfg, grid_config_2d_ptr const &grid_cfg,
                           std::size_t const &y_index, double const &y,
                           boundary_2d_pair const &horizontal_boundary_pair, matrix_2d const &input,
                           matrix_2d const &inhom_input, double const &time, container_t &solution);
};

/**
    heston_boundary_solver object
 */
class heston_boundary_solver
{

  private:
    heat_coefficients_2d_ptr coefficients_;
    grid_config_2d_ptr grid_cfg_;

    explicit heston_boundary_solver() = delete;

  public:
    explicit heston_boundary_solver(heat_coefficients_2d_ptr const &coefficients,
                                    grid_config_2d_ptr const &grid_config);

    ~heston_boundary_solver();

    heston_boundary_solver(heston_boundary_solver const &) = delete;
    heston_boundary_solver(heston_boundary_solver &&) = delete;
    heston_boundary_solver &operator=(heston_boundary_solver const &) = delete;
    heston_boundary_solver &operator=(heston_boundary_solver &&) = delete;

    void solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizonatal_boundary_pair,
               boundary_2d_ptr const &vertical_upper_boundary_ptr, double const &time, matrix_2d &solution);

    void solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizonatal_boundary_pair, double const &time,
               matrix_2d &solution);
};

using heston_boundary_solver_ptr = sptr_t<heston_boundary_solver>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_BOUNDARY_SOLVER_HPP_
