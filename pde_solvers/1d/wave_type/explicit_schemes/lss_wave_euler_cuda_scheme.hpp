/**

    @file      lss_wave_euler_cuda_scheme.hpp
    @brief     Euler CUDA scheme for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_EULER_CUDA_SCHEME_HPP_)
#define _LSS_WAVE_EULER_CUDA_SCHEME_HPP_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../explicit_coefficients/lss_wave_explicit_coefficients.hpp"
#include "../solver_method/lss_wave_euler_cuda_solver_method.hpp"
#include "../time_loop/lss_wave_explicit_time_loop.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_containers::matrix_2d;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;
using lss_utility::NaN;

/**
 * wave_euler_cuda_scheme object
 */

class wave_euler_cuda_scheme
{

  private:
    wave_explicit_coefficients_ptr euler_coeffs_;
    boundary_1d_pair boundary_pair_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    grid_config_1d_ptr grid_cfg_;

    bool is_stable();

    void initialize()
    {
        LSS_ASSERT(is_stable() == true, "The chosen scheme is not stable");
    }

    explicit wave_euler_cuda_scheme() = delete;

  public:
    wave_euler_cuda_scheme(wave_explicit_coefficients_ptr const &coefficients, boundary_1d_pair const &boundary_pair,
                           pde_discretization_config_1d_ptr const &discretization_config,
                           grid_config_1d_ptr const &grid_config);

    ~wave_euler_cuda_scheme();

    void operator()(container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution,
                    bool is_wave_sourse_set, std::function<double(double, double)> const &wave_source,
                    traverse_direction_enum traverse_dir);

    void operator()(container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution,
                    bool is_wave_sourse_set, std::function<double(double, double)> const &wave_source,
                    traverse_direction_enum traverse_dir, matrix_2d &solutions);
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_EULER_CUDA_SCHEME_HPP_
