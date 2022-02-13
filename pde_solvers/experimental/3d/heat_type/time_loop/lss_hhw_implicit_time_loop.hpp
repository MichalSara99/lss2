#if !defined(_LSS_HHW_IMPLICIT_TIME_LOOP_HPP_)
#define _LSS_HHW_IMPLICIT_TIME_LOOP_HPP_

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_enumerations.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../containers/lss_container_2d.hpp"
#include "../../../../../containers/lss_container_3d.hpp"
#include "../../../../../discretization/lss_grid_config.hpp"
#include "../boundary_solver/lss_hhw_explicit_boundary_solver.hpp"
#include "../splitting_method/lss_heat_splitting_method_3d.hpp"

namespace lss_pde_solvers
{

namespace three_dimensional
{

using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_containers::container_3d;
using lss_enumerations::by_enum;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_3d_ptr;
using lss_utility::range_ptr;

/**
    hhw_implicit_boundaries object
 */
class hhw_implicit_boundaries
{

  public:
    static boundary_3d_pair get_y(grid_config_2d_ptr const &grid_config_xz,
                                  container_3d<by_enum::RowPlane> const &next_solution);

    static boundary_3d_pair get_intermed_x(grid_config_2d_ptr const &grid_config_yz,
                                           container_3d<by_enum::RowPlane> const &next_solution);

    static boundary_3d_pair get_intermed_z(grid_config_2d_ptr const &grid_config_xy,
                                           container_3d<by_enum::RowPlane> const &next_solution);
};

/**
 * hhw_implicit_time_loop object
 */
class hhw_implicit_time_loop
{

  public:
    static void run(heat_splitting_method_3d_ptr const &solver_ptr,
                    hhw_explicit_boundary_solver_ptr const &boundary_solver_ptr,
                    boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                    boundary_3d_pair const &z_boundary_pair, grid_config_3d_ptr const &grid_config,
                    range_ptr const &time_range, std::size_t const &last_time_idx, double const time_step,
                    traverse_direction_enum const &traverse_dir, container_3d<by_enum::RowPlane> &prev_solution,
                    container_3d<by_enum::RowPlane> &next_solution);
};
} // namespace three_dimensional
} // namespace lss_pde_solvers

#endif //_LSS_HHW_IMPLICIT_TIME_LOOP_HPP_
