/**

    @file      lss_wave_implicit_time_loop.hpp
    @brief     Implicit time loop for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_IMPLICIT_TIME_LOOP_HPP_)
#define _LSS_WAVE_IMPLICIT_TIME_LOOP_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_range.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../solver_method/lss_wave_implicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_containers::matrix_2d;
using lss_enumerations::traverse_direction_enum;
using lss_utility::container_t;
using lss_utility::range_ptr;

/**
 * wave_implicit_time_loop object
 */
class wave_implicit_time_loop
{

  public:
    static void run(wave_implicit_solver_method_ptr const &solver_ptr, boundary_1d_pair const &boundary_pair,
                    range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                    traverse_direction_enum const &traverse_dir, container_t &prev_solution_0,
                    container_t &prev_solution_1, container_t &next_solution);

    static void run(wave_implicit_solver_method_ptr const &solver_ptr, boundary_1d_pair const &boundary_pair,
                    range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                    traverse_direction_enum const &traverse_dir, container_t &prev_solution_0,
                    container_t &prev_solution_1, std::function<double(double, double)> const &wave_source,
                    container_t &next_solution);

    static void run_with_stepping(wave_implicit_solver_method_ptr const &solver_ptr,
                                  boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                  std::size_t const &last_time_idx, double const time_step,
                                  traverse_direction_enum const &traverse_dir, container_t &prev_solution_0,
                                  container_t &prev_solution_1, container_t &next_solution, matrix_2d &solutions);

    static void run_with_stepping(wave_implicit_solver_method_ptr const &solver_ptr,
                                  boundary_1d_pair const &boundary_pair, range_ptr const &time_range,
                                  std::size_t const &last_time_idx, double const time_step,
                                  traverse_direction_enum const &traverse_dir, container_t &prev_solution_0,
                                  container_t &prev_solution_1,
                                  std::function<double(double, double)> const &wave_source, container_t &next_solution,
                                  matrix_2d &solutions);
};

} // namespace one_dimensional
} // namespace lss_pde_solvers

#endif //_LSS_WAVE_IMPLICIT_TIME_LOOP_HPP_
