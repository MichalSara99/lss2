/**

    @file      lss_ode_equation_implicit_kernel.hpp
    @brief     Implicit kernel for 2nd degree ODEs
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_ODE_EQUATION_IMPLICIT_KERNEL_HPP_)
#define _LSS_ODE_EQUATION_IMPLICIT_KERNEL_HPP_

#include "../../boundaries/lss_boundary.hpp"
#include "../../common/lss_utility.hpp"
#include "../../discretization/lss_discretization.hpp"
#include "../../discretization/lss_grid_config.hpp"
#include "../../ode_solvers/lss_ode_discretization_config.hpp"
#include "../../ode_solvers/lss_ode_solver_config.hpp"
#include "../../ode_solvers/second_degree/implicit_coefficients/lss_ode_implicit_coefficients.hpp"
#include "../../ode_solvers/second_degree/solver_method/lss_ode_implicit_solver_method.hpp"
#include "../../ode_solvers/transformation/lss_ode_data_transform.hpp"
#include "../../sparse_solvers/tridiagonal/lss_tridiagonal_solver.hpp"

namespace lss_ode_solvers
{

using lss_boundary::boundary_1d_pair;
using lss_grids::grid_config_1d_ptr;
using lss_tridiagonal_solver::tridiagonal_solver_ptr;
using lss_utility::container_t;

/**

    @class   ode_equation_implicit_kernel
    @brief   ode_equation_implicit_kernel object
    @details ~

**/
class ode_equation_implicit_kernel
{

  private:
    tridiagonal_solver_ptr solver_;
    boundary_1d_pair boundary_pair_;
    ode_data_transform_ptr ode_data_cfg_;
    ode_discretization_config_ptr discretization_cfg_;
    ode_implicit_solver_config_ptr solver_cfg_;
    grid_config_1d_ptr grid_cfg_;

  public:
    ode_equation_implicit_kernel(tridiagonal_solver_ptr solver, boundary_1d_pair const &boundary_pair,
                                 ode_data_transform_ptr const &ode_data_config,
                                 ode_discretization_config_ptr const &discretization_config,
                                 ode_implicit_solver_config_ptr const &solver_config,
                                 grid_config_1d_ptr const &grid_config);

    void operator()(container_t &solution, bool is_ode_nonhom_set, std::function<double(double)> const &ode_nonhom,
                    double omega_value);
};

} // namespace lss_ode_solvers

#endif ///_LSS_ODE_EQUATION_IMPLICIT_KERNEL_HPP_
