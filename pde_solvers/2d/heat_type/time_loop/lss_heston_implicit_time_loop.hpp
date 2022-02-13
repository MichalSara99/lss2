#if !defined(_LSS_HESTON_IMPLICIT_TIME_LOOP_HPP_)
#define _LSS_HESTON_IMPLICIT_TIME_LOOP_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_container_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../boundary_solver/lss_heston_explicit_boundary_solver.hpp"
#include "../explicit_coefficients/lss_heston_euler_coefficients.hpp"
#include "../splitting_method/lss_heat_splitting_method.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::container_2d;
using lss_enumerations::by_enum;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::range_ptr;

/**
    heston_implicit_boundaries object
 */
class heston_implicit_boundaries
{

  public:
    static boundary_2d_pair get_vertical(grid_config_1d_ptr const &grid_config_x,
                                         container_2d<by_enum::Row> const &next_solution);

    static boundary_2d_pair get_intermed_horizontal(grid_config_1d_ptr const &grid_config_y,
                                                    container_2d<by_enum::Row> const &next_solution);
};

/**
 * heston_implicit_time_loop object
 */
class heston_implicit_time_loop
{

  public:
    static void run(heat_splitting_method_ptr const &solver_ptr,
                    heston_explicit_boundary_solver_ptr const &boundary_solver_ptr,
                    boundary_2d_pair const &horizontal_boundary_pair,
                    boundary_2d_ptr const &vertical_upper_boundary_ptr, grid_config_2d_ptr const &grid_config,
                    range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                    traverse_direction_enum const &traverse_dir, container_2d<by_enum::Row> &prev_solution,
                    container_2d<by_enum::Row> &next_solution);
};
} // namespace two_dimensional
} // namespace lss_pde_solvers

#endif //_LSS_HESTON_IMPLICIT_TIME_LOOP_HPP_
