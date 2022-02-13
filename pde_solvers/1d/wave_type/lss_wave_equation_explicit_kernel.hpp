/**

    @file      lss_wave_equation_explicit_kernel.hpp
    @brief     Explicit kernel for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_EQUATION_EXPLICIT_KERNEL_HPP_)
#define _LSS_WAVE_EQUATION_EXPLICIT_KERNEL_HPP_

#include <vector>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../lss_wave_data_config.hpp"
#include "../../lss_wave_solver_config.hpp"
#include "../../transformation/lss_wave_data_transform.hpp"
#include "explicit_coefficients/lss_wave_explicit_coefficients.hpp"
#include "explicit_schemes/lss_wave_euler_cuda_scheme.hpp"
#include "explicit_schemes/lss_wave_euler_scheme.hpp"

namespace lss_pde_solvers
{
namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_containers::matrix_2d;
using lss_enumerations::dimension_enum;
using lss_enumerations::explicit_pde_schemes_enum;
using lss_enumerations::memory_space_enum;
using lss_utility::container_t;

template <memory_space_enum memory_enum> class wave_equation_explicit_kernel
{
};

// ===================================================================
// ============================== DEVICE =============================
// ===================================================================
template <> class wave_equation_explicit_kernel<memory_space_enum::Device>
{

  private:
    boundary_1d_pair boundary_pair_;
    wave_data_transform_1d_ptr wave_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    wave_explicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    wave_equation_explicit_kernel(boundary_1d_pair const &boundary_pair,
                                  wave_data_transform_1d_ptr const &wave_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  wave_explicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution,
                    bool is_wave_sourse_set, std::function<double(double, double)> const &wave_source);

    void operator()(container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution,
                    bool is_wave_sourse_set, std::function<double(double, double)> const &wave_source,
                    matrix_2d &solutions);
};

// ===================================================================
// ============================== HOST ===============================
// ===================================================================

template <> class wave_equation_explicit_kernel<memory_space_enum::Host>
{

  private:
    boundary_1d_pair boundary_pair_;
    wave_data_transform_1d_ptr wave_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    wave_explicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    wave_equation_explicit_kernel(boundary_1d_pair const &boundary_pair,
                                  wave_data_transform_1d_ptr const &wave_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  wave_explicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution,
                    bool is_wave_sourse_set, std::function<double(double, double)> const &wave_source);

    void operator()(container_t &prev_solution_0, container_t &prev_solution_1, container_t &next_solution,
                    bool is_wave_sourse_set, std::function<double(double, double)> const &wave_source,
                    matrix_2d &solutions);
};

} // namespace one_dimensional
} // namespace lss_pde_solvers
#endif ///_LSS_WAVE_EQUATION_EXPLICIT_KERNEL_HPP_
