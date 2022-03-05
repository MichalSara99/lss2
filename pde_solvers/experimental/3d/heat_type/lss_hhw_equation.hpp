/**

    @file      lss_hhw_equation.hpp
    @brief     Heston-Hull-White equation type
    @details   ~
    @author    Michal Sara
    @date      5.03.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HHW_EQUATION_HPP_)
#define _LSS_HHW_EQUATION_HPP_

#include <functional>
#include <map>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../containers/lss_matrix_3d.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config_hints.hpp"
#include "../../../../discretization/lss_grid_transform_config.hpp"
#include "../../../lss_heat_data_config.hpp"
#include "../../../lss_heat_solver_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../../../lss_splitting_method_config.hpp"
#include "../../../transformation/lss_heat_data_transform.hpp"
#include "transformation/lss_hhw_boundary_transform.hpp"

namespace lss_pde_solvers
{

namespace three_dimensional
{
using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_containers::matrix_3d;
using lss_enumerations::grid_enum;
using lss_grids::grid_config_hints_3d_ptr;
using lss_grids::grid_transform_config_3d_ptr;

using d_3d = discretization_3d;

namespace implicit_solvers
{

/**

    @class   hhw_equation
    @brief   Represents general variable coefficient Heston-Hull-White type equation
    @details

            u_t = a(t,x,y,z)*u_xx + b(t,x,y,z)*u_yy + c(t,x,y,z)*u_zz
                + d(t,x,y,z)*u_xy + e(t,x,y,z)*u_xz + f(t,x,y,z)*u_yz
                + g(t,x,y,z)*u_x + h(t,x,y,z)*u_y + i(t,x,y,z)*u_z
                + j(t,x,y,z)*u + R(t,x,y,z)

            t > 0, x_1 < x < x_2, y_1 < y < y_2, z_1 < z < z_2

            with initial condition:

            u(0,x,y,z) = G(x,y,z)

            or terminal condition:

            u(T,x,y,z) = G(x,y,z)

            horizontal_boundary_pair = S = (S_1,S_2) boundary

                       _________________
                      /                /|
                     /    vol (Y)     / |
                    /________________/  |
                    |S_1,S_1,S_1,S_1|   |
                    |               |   |
                    |               |   |
            S (X)   |               |   |
                    |               |   |
                    |               |   |
                    |               |   /
                    |               |  / rate (Z)
                    |S_2,S_2,S_2,S_2| /
                    |_______________|/



**/
class hhw_equation
{

  private:
    heat_data_transform_3d_ptr heat_data_trans_cfg_;
    pde_discretization_config_3d_ptr discretization_cfg_;
    hhw_boundary_transform_ptr hhw_boundary_;
    splitting_method_config_ptr splitting_method_cfg_;
    grid_transform_config_3d_ptr grid_trans_cfg_; // this may be removed as it is not used later
    heat_implicit_solver_config_ptr solver_cfg_;
    std::map<std::string, double> solver_config_details_;

    explicit hhw_equation() = delete;

    void initialize(heat_data_config_3d_ptr const &heat_data_cfg, grid_config_hints_3d_ptr const &grid_config_hints,
                    boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                    boundary_3d_pair const &z_boundary_pair);

  public:
    explicit hhw_equation(heat_data_config_3d_ptr const &heat_data_config,
                          pde_discretization_config_3d_ptr const &discretization_config,
                          boundary_3d_pair const &x_boundary_pair, boundary_3d_ptr const &y_upper_boundary_ptr,
                          boundary_3d_pair const &z_boundary_pair,
                          splitting_method_config_ptr const &splitting_method_config,
                          grid_config_hints_3d_ptr const &grid_config_hints,
                          heat_implicit_solver_config_ptr const &solver_config =
                              default_heat_solver_configs::host_fwd_dssolver_euler_solver_config_ptr,
                          std::map<std::string, double> const &solver_config_details = std::map<std::string, double>());

    ~hhw_equation();

    hhw_equation(hhw_equation const &) = delete;
    hhw_equation(hhw_equation &&) = delete;
    hhw_equation &operator=(hhw_equation const &) = delete;
    hhw_equation &operator=(hhw_equation &&) = delete;

    /**
     * Get the final solution of the PDE
     *
     * \param solution -  3D container for solution
     */
    LSS_API void solve(matrix_3d &solution);
};

} // namespace implicit_solvers

} // namespace three_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HHW_EQUATION_HPP_
