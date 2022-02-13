/**

    @file      lss_heat_equation_implicit_kernel.hpp
    @brief     Implicit heat kernel
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EQUATION_IMPLICIT_KERNEL_HPP_)
#define _LSS_HEAT_EQUATION_IMPLICIT_KERNEL_HPP_

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../lss_heat_solver_config.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../transformation/lss_heat_data_transform.hpp"
#include "implicit_coefficients/lss_heat_coefficients.hpp"

namespace lss_pde_solvers
{
namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_containers::matrix_2d;
using lss_enumerations::dimension_enum;
using lss_enumerations::implicit_pde_schemes_enum;
using lss_enumerations::memory_space_enum;
using lss_enumerations::traverse_direction_enum;
using lss_enumerations::tridiagonal_method_enum;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;

template <memory_space_enum memory_enum, tridiagonal_method_enum tridiagonal_method> class heat_equation_implicit_kernel
{
};

// ===================================================================
// ============================== DEVICE =============================
// ===================================================================
template <> class heat_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::CUDASolver>
{

  private:
    boundary_1d_pair boundary_pair_;
    heat_data_transform_1d_ptr heat_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    heat_equation_implicit_kernel(boundary_1d_pair const &boundary_pair,
                                  heat_data_transform_1d_ptr const &heat_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  heat_implicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, matrix_2d &solutions);
};

template <> class heat_equation_implicit_kernel<memory_space_enum::Device, tridiagonal_method_enum::SORSolver>
{

  private:
    boundary_1d_pair boundary_pair_;
    heat_data_transform_1d_ptr heat_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    heat_equation_implicit_kernel(boundary_1d_pair const &boundary_pair,
                                  heat_data_transform_1d_ptr const &heat_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  heat_implicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, double omega_value);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, double omega_value, matrix_2d &solutions);
};

// ===================================================================
// ================================ HOST =============================
// ===================================================================
template <> class heat_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::CUDASolver>
{

  private:
    boundary_1d_pair boundary_pair_;
    heat_data_transform_1d_ptr heat_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    heat_equation_implicit_kernel(boundary_1d_pair const &boundary_pair,
                                  heat_data_transform_1d_ptr const &heat_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  heat_implicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, matrix_2d &solutions);
};

template <> class heat_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::SORSolver>
{
  private:
    boundary_1d_pair boundary_pair_;
    heat_data_transform_1d_ptr heat_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    heat_equation_implicit_kernel(boundary_1d_pair const &boundary_pair,
                                  heat_data_transform_1d_ptr const &heat_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  heat_implicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, double omega_value);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, double omega_value, matrix_2d &solutions);
};

template <> class heat_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::DoubleSweepSolver>
{

  private:
    boundary_1d_pair boundary_pair_;
    heat_data_transform_1d_ptr heat_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    heat_equation_implicit_kernel(boundary_1d_pair const &boundary_pair,
                                  heat_data_transform_1d_ptr const &heat_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  heat_implicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, matrix_2d &solutions);
};

template <> class heat_equation_implicit_kernel<memory_space_enum::Host, tridiagonal_method_enum::ThomasLUSolver>
{

  private:
    boundary_1d_pair boundary_pair_;
    heat_data_transform_1d_ptr heat_data_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    heat_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    heat_equation_implicit_kernel(boundary_1d_pair const &boundary_pair,
                                  heat_data_transform_1d_ptr const &heat_data_config,
                                  pde_discretization_config_1d_ptr const &discretization_config,
                                  heat_implicit_solver_config_ptr const &solver_config,
                                  grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source);

    void operator()(container_t &solution, bool is_heat_sourse_set,
                    std::function<double(double, double)> const &heat_source, matrix_2d &solutions);
};
} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EQUATION_IMPLICIT_KERNEL_HPP_
