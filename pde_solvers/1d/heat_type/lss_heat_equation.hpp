/**

    @file      lss_heat_equation.hpp
    @brief     Heat equation solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_1D_GENERAL_HEAT_EQUATION_HPP_)
#define _LSS_1D_GENERAL_HEAT_EQUATION_HPP_

#include <functional>
#include <map>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../discretization/lss_grid.hpp"
#include "../../../discretization/lss_grid_config.hpp"
#include "../../../discretization/lss_grid_config_hints.hpp"
#include "../../../discretization/lss_grid_transform_config.hpp"
#include "../../../transformation/lss_boundary_transform.hpp"
#include "../../lss_heat_data_config.hpp"
#include "../../lss_heat_solver_config.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../transformation/lss_heat_data_transform.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{
using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_containers::matrix_2d;
using lss_enumerations::grid_enum;
using lss_grids::grid_config_1d;
using lss_grids::grid_config_hints_1d_ptr;
using lss_grids::grid_transform_config_1d;
using lss_grids::grid_transform_config_1d_ptr;
using lss_transformation::boundary_transform_1d;
using lss_transformation::boundary_transform_1d_ptr;
using lss_utility::container_t;

namespace implicit_solvers
{

/**
    @class   heat_equation
    @brief   Represents general variable coefficient 1D heat equation solver
    @details Heat equation is of following standard form:

             u_t = a(t,x)*u_xx + b(t,x)*u_x + c(t,x)*u + F(t,x),
             x_1 < x < x_2
             t_1 < t < t_2

             with initial condition:

                u(t_1,x) = f(x)

             or terminal condition:

                u(t_2,x) = f(x)

**/
class heat_equation
{

  private:
    heat_data_transform_1d_ptr heat_data_trans_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    boundary_transform_1d_ptr boundary_;
    grid_transform_config_1d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    heat_implicit_solver_config_ptr solver_cfg_;
    std::map<std::string, double> solver_config_details_;

    explicit heat_equation() = delete;

    void initialize(heat_data_config_1d_ptr const &heat_data_cfg, grid_config_hints_1d_ptr const &grid_config_hints,
                    boundary_1d_pair const &boundary_pair);

  public:
    explicit heat_equation(
        heat_data_config_1d_ptr const &heat_data_config, pde_discretization_config_1d_ptr const &discretization_config,
        boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
        heat_implicit_solver_config_ptr const &solver_config =
            default_heat_solver_configs::host_fwd_dssolver_euler_solver_config_ptr,
        std::map<std::string, double> const &solver_config_details = std::map<std::string, double>());

    ~heat_equation();

    heat_equation(heat_equation const &) = delete;
    heat_equation(heat_equation &&) = delete;
    heat_equation &operator=(heat_equation const &) = delete;
    heat_equation &operator=(heat_equation &&) = delete;

    /**
        @brief  Get the final solution of the PDE
        @param solution - container for solution
        @retval
    **/
    LSS_API void solve(container_t &solution);

    /**
        @brief  Get all solutions in time (surface) of the PDE
        @param solutions - matrix_2d containing all solutions (one solution per row)
        @retval
    **/
    LSS_API void solve(matrix_2d &solutions);
};

} // namespace implicit_solvers

namespace explicit_solvers
{

/**

    @class   heat_equation
    @brief   Represents general variable coefficient 1D heat equation solver
    @details Heat equation is of following standard form:

             u_t = a(t,x)*u_xx + b(t,x)*u_x + c(t,x)*u + F(t,x),
             x_1 < x < x_2
             t_1 < t < t_2

             with initial condition:

                u(t_1,x) = f(x)

             or terminal condition:

                u(t_2,x) = f(x)

**/
class heat_equation
{
  private:
    heat_data_transform_1d_ptr heat_data_trans_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    boundary_transform_1d_ptr boundary_;
    grid_transform_config_1d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    heat_explicit_solver_config_ptr solver_cfg_;

    explicit heat_equation() = delete;

    void initialize(heat_data_config_1d_ptr const &heat_data_cfg, grid_config_hints_1d_ptr const &grid_config_hints,
                    boundary_1d_pair const &boundary_pair);

  public:
    explicit heat_equation(heat_data_config_1d_ptr const &heat_data_config,
                           pde_discretization_config_1d_ptr const &discretization_config,
                           boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
                           heat_explicit_solver_config_ptr const &solver_config =
                               default_heat_solver_configs::dev_expl_fwd_euler_solver_config_ptr);

    ~heat_equation();

    heat_equation(heat_equation const &) = delete;
    heat_equation(heat_equation &&) = delete;
    heat_equation &operator=(heat_equation const &) = delete;
    heat_equation &operator=(heat_equation &&) = delete;

    /**
        @brief  Get the final solution of the PDE
        @param  container for solution
        @retval
    **/
    LSS_API void solve(container_t &solution);

    /**
        @brief  Get all solutions in time (surface) of the PDE
        @param  matrix_2d containing all solutions (one solution per row)
        @retval
    **/
    LSS_API void solve(matrix_2d &solutions);
};

} // namespace explicit_solvers

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_1D_GENERAL_HEAT_EQUATION_HPP_
