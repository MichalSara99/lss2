/**

    @file      lss_wave_equation.hpp
    @brief     Wave equation
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_EQUATION_HPP_)
#define _LSS_WAVE_EQUATION_HPP_

#include <functional>
#include <map>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../containers/lss_matrix_2d.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../transformation/lss_boundary_transform.hpp"
#include "../../lss_pde_discretization_config.hpp"
#include "../../lss_wave_data_config.hpp"
#include "../../lss_wave_solver_config.hpp"
#include "../../transformation/lss_wave_data_transform.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{
using lss_boundary::boundary_1d_pair;
using lss_containers::matrix_2d;
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

    @class   wave_equation
    @brief   Represents general variable coefficient 1D wave equation solver
    @details Wave equation is of following standard form:

            u_tt + a(t,x)*u_t = b(t,x)*u_xx + c(t,x)*u_x + d(t,x)*u + F(t,x),
            x_1 < x < x_2
            t_1 < t < t_2

            with initial condition:

            u(t_1,x) = f(x)

            u_t(t_1,x) = g(x)

            or terminal condition:

            u(t_2,x) = f(x)

            u_t(t_2,x) = g(x)

**/
class wave_equation
{

  private:
    wave_data_transform_1d_ptr wave_data_trans_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    boundary_transform_1d_ptr boundary_;
    grid_transform_config_1d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    wave_implicit_solver_config_ptr solver_cfg_;
    std::map<std::string, double> solver_config_details_;

    explicit wave_equation() = delete;

    void initialize(wave_data_config_1d_ptr const &wave_data_cfg, grid_config_hints_1d_ptr const &grid_config_hints,
                    boundary_1d_pair const &boundary_pair);

  public:
    explicit wave_equation(
        wave_data_config_1d_ptr const &wave_data_config, pde_discretization_config_1d_ptr const &discretization_config,
        boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
        wave_implicit_solver_config_ptr const &solver_config =
            default_wave_solver_configs::dev_fwd_cusolver_qr_solver_config_ptr,
        std::map<std::string, double> const &solver_config_details = std::map<std::string, double>());

    ~wave_equation();

    wave_equation(wave_equation const &) = delete;
    wave_equation(wave_equation &&) = delete;
    wave_equation &operator=(wave_equation const &) = delete;
    wave_equation &operator=(wave_equation &&) = delete;

    /**
        @brief  Get the final solution of the PDE
        @param  solution - container for solution
        @retval
    **/
    LSS_API void solve(container_t &solution);

    /**
        @brief  Get all solutions in time (surface) of the PDE
        @param  solution - matrix_2d containing all solutions (one solution per row)
        @retval
    **/
    LSS_API void solve(matrix_2d &solutions);
};

} // namespace implicit_solvers

namespace explicit_solvers
{

/**

    @class   wave_equation
    @brief   Represents general variable coefficient 1D wave equation solver
    @details Wave equation is of following standard form:

            u_tt + a(t,x)*u_t = b(t,x)*u_xx + c(t,x)*u_x + d(t,x)*u + F(t,x),
            x_1 < x < x_2
            t_1 < t < t_2

            with initial condition:

            u(t_1,x) = f(x)

            u_t(t_1,x) = g(x)

            or terminal condition:

            u(t_2,x) = f(x)

            u_t(t_2,x) = g(x)

**/
class wave_equation
{
  private:
    wave_data_transform_1d_ptr wave_data_trans_cfg_;
    pde_discretization_config_1d_ptr discretization_cfg_;
    boundary_transform_1d_ptr boundary_;
    grid_transform_config_1d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    wave_explicit_solver_config_ptr solver_cfg_;

    explicit wave_equation() = delete;

    void initialize(wave_data_config_1d_ptr const &wave_data_cfg, grid_config_hints_1d_ptr const &grid_config_hints,
                    boundary_1d_pair const &boundary_pair);

  public:
    explicit wave_equation(wave_data_config_1d_ptr const &wave_data_config,
                           pde_discretization_config_1d_ptr const &discretization_config,
                           boundary_1d_pair const &boundary_pair, grid_config_hints_1d_ptr const &grid_config_hints,
                           wave_explicit_solver_config_ptr const &solver_config =
                               default_wave_solver_configs::dev_expl_fwd_solver_config_ptr);

    ~wave_equation();

    wave_equation(wave_equation const &) = delete;
    wave_equation(wave_equation &&) = delete;
    wave_equation &operator=(wave_equation const &) = delete;
    wave_equation &operator=(wave_equation &&) = delete;

    /**
        @brief  Get the final solution of the PDE
        @param  solution - container for solution
        @retval
    **/
    LSS_API void solve(container_t &solution);

    /**
        @brief  Get all solutions in time (surface) of the PDE
        @param  solution - matrix_2d containing all solutions (one solution per row)
        @retval
    **/
    LSS_API void solve(matrix_2d &solutions);
};

} // namespace explicit_solvers

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_EQUATION_HPP_
