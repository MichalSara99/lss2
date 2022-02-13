/**

    @file      lss_ode_equation.hpp
    @brief     2nd degree ODE equation solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_ODE_EQUATION_HPP_)
#define _LSS_ODE_EQUATION_HPP_

#include <map>

#include "../../boundaries/lss_boundary.hpp"
#include "../../common/lss_macros.hpp"
//#include "lss_general_ode_equation_explicit_kernel.hpp"
#include "../../discretization/lss_grid_config.hpp"
#include "../../ode_solvers/lss_ode_data_config.hpp"
#include "../../ode_solvers/lss_ode_discretization_config.hpp"
#include "../../ode_solvers/lss_ode_solver_config.hpp"
#include "../../ode_solvers/transformation/lss_ode_data_transform.hpp"
#include "../../transformation/lss_boundary_transform.hpp"
#include "lss_ode_equation_implicit_kernel.hpp"

namespace lss_ode_solvers
{

using lss_boundary::boundary_1d_pair;
using lss_grids::grid_config_hints_1d_ptr;
using lss_grids::grid_transform_config_1d_ptr;
using lss_transformation::boundary_transform_1d_ptr;
using lss_utility::container_t;

namespace implicit_solvers
{
/**

    @class   ode_equation
    @brief   Represents general 2.degree ODE
    @details General 2.degree ODE of the following form

    u''(x) + a(x)*u'(x) + b(x)*u(x) = g(x),
    x_1 < x < x_2

**/
class ode_equation
{

  private:
    ode_data_transform_ptr ode_data_trans_cfg_;
    ode_discretization_config_ptr ode_discretization_cfg_;
    boundary_transform_1d_ptr boundary_;
    grid_transform_config_1d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    ode_implicit_solver_config_ptr ode_solver_cfg_;
    std::map<std::string, double> ode_solver_config_details_;

    explicit ode_equation() = delete;

    void initialize(ode_data_config_ptr const &ode_data_cfg, grid_config_hints_1d_ptr const &grid_config_hints,
                    boundary_1d_pair const &boundary_pair);

  public:
    explicit ode_equation(
        ode_data_config_ptr const &ode_data_config, ode_discretization_config_ptr const &ode_discretization_config,
        boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
        ode_implicit_solver_config_ptr const &ode_solver_config =
            default_ode_solver_configs::dev_cusolver_qr_solver_config_ptr,
        std::map<std::string, double> const &ode_solver_config_details = std::map<std::string, double>());

    ~ode_equation();

    ode_equation(ode_equation const &) = delete;
    ode_equation(ode_equation &&) = delete;
    ode_equation &operator=(ode_equation const &) = delete;
    ode_equation &operator=(ode_equation &&) = delete;

    /**
        @brief Get the final solution of the ODE
        @param solution - container for solution
    **/
    LSS_API void solve(container_t &solution);
};

} // namespace implicit_solvers

namespace explicit_solvers
{

} // namespace explicit_solvers

} // namespace lss_ode_solvers

#endif ///_LSS_ODE_EQUATION_HPP_
