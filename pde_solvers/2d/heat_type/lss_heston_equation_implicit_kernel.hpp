/**

    @file      lss_heston_equation_implicit_kernel.hpp
    @brief     Implicit kernel for Heston type problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright ? Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HESTON_EQUATION_IMPLICIT_KERNEL_HPP_)
#define _LSS_HESTON_EQUATION_IMPLICIT_KERNEL_HPP_

#include <vector>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../containers/lss_matrix_3d.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../discretization/lss_grid.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../lss_heat_solver_config.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../lss_splitting_method_config.hpp"
#include "../../transformation/lss_heat_data_transform.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_boundary::boundary_2d_pair;
using lss_boundary::boundary_2d_ptr;
using lss_containers::matrix_2d;
using lss_containers::matrix_3d;
using lss_enumerations::memory_space_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_grids::grid_config_2d_ptr;
using lss_utility::sptr_t;

template <memory_space_enum memory_enum, tridiagonal_method_enum tridiagonal_method>
class heston_equation_implicit_kernel
{
};

// ===================================================================
// ============================== DEVICE =============================
// ===================================================================
template <> class heston_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>
{

  private:
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    heat_data_transform_2d_ptr heat_data_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    splitting_method_config_ptr splitting_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_2d_ptr grid_cfg_;

  public:
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, matrix_3d &solutions);
};

template <> class heston_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>
{

  private:
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    heat_data_transform_2d_ptr heat_data_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    splitting_method_config_ptr splitting_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_2d_ptr grid_cfg_;

  public:
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, double omega_value);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, double omega_value,
                    matrix_3d &solutions);
};

// ===================================================================
// ================================ HOST =============================
// ===================================================================
template <> class heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>
{

  private:
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    heat_data_transform_2d_ptr heat_data_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    splitting_method_config_ptr splitting_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_2d_ptr grid_cfg_;

  public:
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, matrix_3d &solutions);
};

template <> class heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>
{
  private:
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    heat_data_transform_2d_ptr heat_data_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    splitting_method_config_ptr splitting_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_2d_ptr grid_cfg_;

  public:
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, double omega_value);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, double omega_value,
                    matrix_3d &solutions);
};

template <> class heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>
{

  private:
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    heat_data_transform_2d_ptr heat_data_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    splitting_method_config_ptr splitting_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_2d_ptr grid_cfg_;

  public:
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, matrix_3d &solutions);
};

template <> class heston_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>
{

  private:
    boundary_2d_ptr boundary_ver_;
    boundary_2d_pair boundary_pair_hor_;
    heat_data_transform_2d_ptr heat_data_cfg_;
    pde_discretization_config_2d_ptr discretization_cfg_;
    splitting_method_config_ptr splitting_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_2d_ptr grid_cfg_;

  public:
    heston_equation_implicit_kernel(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                    boundary_2d_pair const &horizontal_boundary_pair,
                                    heat_data_transform_2d_ptr const &heat_data_config,
                                    pde_discretization_config_2d_ptr const &discretization_config,
                                    splitting_method_config_ptr const &splitting_config,
                                    heat_implicit_solver_config_ptr const &solver_config,
                                    grid_config_2d_ptr const &grid_config);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source);

    void operator()(matrix_2d &prev_solution, matrix_2d &next_solution, bool is_heat_sourse_set,
                    std::function<double(double, double, double)> const &heat_source, matrix_3d &solutions);
};
} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EQUATION_IMPLICIT_KERNEL_HPP_
