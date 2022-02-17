/**

    @file      lss_implicit_time_loop_2d.hpp
    @brief     Implicit time loop for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_IMPLICIT_TIME_LOOP_2D_HPP_)
#define _LSS_IMPLICIT_TIME_LOOP_2D_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../boundary_solver/lss_heston_boundary_solver.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients_2d.hpp"
#include "../solver_method/splitting/lss_heat_splitting_method.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::matrix_2d;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::range_ptr;

/**
    implicit_boundary object
 */
class implicit_boundary
{

  public:
    static boundary_2d_pair get_vertical(grid_config_1d_ptr const &grid_config_x, matrix_2d const &next_solution);

    static boundary_2d_pair get_intermed_horizontal(grid_config_1d_ptr const &grid_config_y,
                                                    matrix_2d const &next_solution);
};

/**
 * implicit_time_loop_2d object
 */
class implicit_time_loop_2d
{

  public:
    static void run(heat_splitting_method_ptr const &solver_ptr, heston_boundary_solver_ptr const &boundary_solver_ptr,
                    boundary_2d_pair const &horizontal_boundary_pair,
                    boundary_2d_ptr const &vertical_upper_boundary_ptr, grid_config_2d_ptr const &grid_config,
                    range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                    traverse_direction_enum const &traverse_dir, matrix_2d &prev_solution, matrix_2d &next_solution);
};
} // namespace two_dimensional
} // namespace lss_pde_solvers

#endif //_LSS_IMPLICIT_TIME_LOOP_2D_HPP_
