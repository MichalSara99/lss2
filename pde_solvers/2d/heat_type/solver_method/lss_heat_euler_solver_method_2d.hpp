/**

    @file      lss_heat_euler_solver_method_2d.hpp
    @brief     Explicit Euler method for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_SOLVER_METHOD_2D_HPP_)
#define _LSS_HEAT_EULER_SOLVER_METHOD_2D_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients_2d.hpp"
#include "lss_heat_explicit_solver_method_2d.hpp"

namespace lss_pde_solvers
{

using lss_containers::matrix_2d;
using lss_containers::matrix_2d_ptr;
using lss_grids::grid_config_2d_ptr;

namespace two_dimensional
{

class explicit_euler_rhs
{

  public:
    static void rhs(heat_euler_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                    matrix_2d const &input, double const &time, matrix_2d &solution);

    static void rhs_source(heat_euler_coefficients_2d_ptr const &cfs, grid_config_2d_ptr const &grid_config,
                           matrix_2d const &input, double const &time, matrix_2d const &inhom_input,
                           matrix_2d &solution);
};

/**
    heat_euler_solver_method_2d object
*/

class heat_euler_solver_method_2d : public heat_explicit_solver_method_2d
{

  private:
    // scheme coefficients:
    heat_euler_coefficients_2d_ptr coefficients_;
    matrix_2d_ptr source_;

    explicit heat_euler_solver_method_2d() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heat_euler_solver_method_2d(heat_euler_coefficients_2d_ptr const &coefficients,
                                         grid_config_2d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_euler_solver_method_2d();

    heat_euler_solver_method_2d(heat_euler_solver_method_2d const &) = delete;
    heat_euler_solver_method_2d(heat_euler_solver_method_2d &&) = delete;
    heat_euler_solver_method_2d &operator=(heat_euler_solver_method_2d const &) = delete;
    heat_euler_solver_method_2d &operator=(heat_euler_solver_method_2d &&) = delete;

    void solve(matrix_2d &prev_solution, double const &time, matrix_2d &solution) override;

    void solve(matrix_2d &prev_solution, double const &time,
               std::function<double(double, double, double)> const &heat_source, matrix_2d &solution) override;
};

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EULER_SOLVER_METHOD_2D_HPP_
