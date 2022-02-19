/**

    @file      lss_heat_euler_cuda_scheme_2d.hpp
    @brief     Euler CUDA scheme for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_CUDA_SCHEME_2D_HPP_)
#define _LSS_HEAT_EULER_CUDA_SCHEME_2D_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../containers/lss_matrix_3d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../boundary_solver/lss_heston_boundary_solver.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients_2d.hpp"
#include "../solver_method/lss_heat_euler_cuda_solver_method_2d.hpp"
#include "../time_loop/lss_explicit_time_loop_2d.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::matrix_2d;
using lss_containers::matrix_3d;
using lss_enumerations::traverse_direction_enum;

class heat_euler_cuda_scheme_2d
{

  private:
    heat_euler_coefficients_2d_ptr euler_coeffs_;
    heston_boundary_solver_ptr heston_boundary_;
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    grid_config_2d_ptr grid_cfg_;

    bool is_stable(heat_coefficients_2d_ptr const &coefficients);

    void initialize(heat_coefficients_2d_ptr const &coefficients);

    explicit heat_euler_cuda_scheme_2d() = delete;

  public:
    heat_euler_cuda_scheme_2d(heat_coefficients_2d_ptr const &coefficients,
                              boundary_2d_ptr const &vertical_upper_boundary_ptr,
                              boundary_2d_pair const &horizontal_boundary_pair,
                              pde_discretization_config_2d_ptr const &discretization_config,
                              grid_config_2d_ptr const &grid_config);

    ~heat_euler_cuda_scheme_2d();

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source,
                    traverse_direction_enum traverse_dir);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source,
                    traverse_direction_enum traverse_dir, matrix_3d &solutions);
};

} // namespace two_dimensional
} // namespace lss_pde_solvers
#endif ///_LSS_HEAT_EULER_CUDA_SCHEME_2D_HPP_
