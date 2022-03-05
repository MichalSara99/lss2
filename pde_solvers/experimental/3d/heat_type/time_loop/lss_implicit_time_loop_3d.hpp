/**

    @file      lss_implicit_time_loop_3d.hpp
    @brief     Implicit time loop for 3D heat problems
    @details   ~
    @author    Michal Sara
    @date      5.03.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_IMPLICIT_TIME_LOOP_3D_HPP_)
#define _LSS_IMPLICIT_TIME_LOOP_3D_HPP_

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_enumerations.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../containers/lss_matrix_2d.hpp"
#include "../../../../../containers/lss_matrix_3d.hpp"
#include "../../../../../discretization/lss_grid_config.hpp"
#include "../boundary_solver/lss_hhw_explicit_boundary_solver.hpp"
#include "../solver_method/splitting/lss_heat_splitting_method_3d.hpp"

namespace lss_pde_solvers
{

namespace three_dimensional
{

using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_containers::matrix_3d;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_2d_ptr;
using lss_grids::grid_config_3d_ptr;
using lss_utility::range_ptr;

/**
    implicit_boundary object
 */
class implicit_boundary
{

  public:
    static boundary_3d_pair get_y(grid_config_2d_ptr const &grid_config_xz, matrix_3d const &next_solution);

    static boundary_3d_pair get_intermed_x(grid_config_2d_ptr const &grid_config_yz, matrix_3d const &next_solution);

    static boundary_3d_pair get_intermed_z(grid_config_2d_ptr const &grid_config_xy, matrix_3d const &next_solution);
};

/**
 * implicit_time_loop_3d object
 */
class implicit_time_loop_3d
{

  public:
    static void run(heat_splitting_method_3d_ptr const &solver_ptr,
                    hhw_explicit_boundary_solver_ptr const &boundary_solver_ptr,
                    boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                    boundary_3d_pair const &z_boundary_pair, grid_config_3d_ptr const &grid_config,
                    range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                    traverse_direction_enum const &traverse_dir, matrix_3d &prev_solution, matrix_3d &next_solution);
};
} // namespace three_dimensional
} // namespace lss_pde_solvers

#endif //_LSS_IMPLICIT_TIME_LOOP_3D_HPP_
