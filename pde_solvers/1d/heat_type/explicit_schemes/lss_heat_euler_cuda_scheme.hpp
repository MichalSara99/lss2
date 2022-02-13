/**

    @file      lss_heat_euler_cuda_scheme.hpp
    @brief     Euler CUDA scheme for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_CUDA_SCHEME_HPP_)
#define _LSS_HEAT_EULER_CUDA_SCHEME_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients.hpp"
#include "../solver_method/lss_heat_euler_cuda_solver_method.hpp"
#include "../time_loop/lss_explicit_time_loop.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_containers::matrix_2d;
using lss_enumerations::traverse_direction_enum;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;
using lss_utility::sptr_t;

/**
 * heat_euler_cuda_scheme object
 */

class heat_euler_cuda_scheme
{

  private:
    heat_euler_coefficients_ptr euler_coeffs_;
    boundary_1d_pair boundary_pair_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    grid_config_1d_ptr grid_cfg_;

    bool is_stable(heat_coefficients_ptr const &coefficients);

    void initialize(heat_coefficients_ptr const &coefficients)
    {
        LSS_ASSERT(is_stable(coefficients) == true, "The chosen scheme is not stable");
        euler_coeffs_ = std::make_shared<heat_euler_coefficients>(coefficients);
    }

    explicit heat_euler_cuda_scheme() = delete;

  public:
    heat_euler_cuda_scheme(heat_coefficients_ptr const &coefficients, boundary_1d_pair const &boundary_pair,
                           pde_discretization_config_1d_ptr const &discretization_config,
                           grid_config_1d_ptr const &grid_config);

    ~heat_euler_cuda_scheme();

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, traverse_direction_enum traverse_dir);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, traverse_direction_enum traverse_dir,
                    matrix_2d &solutions);
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EULER_CUDA_SCHEME_HPP_
