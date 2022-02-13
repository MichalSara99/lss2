/**

    @file      lss_heat_euler_solver_method.hpp
    @brief     Euler solver method for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_SOLVER_METHOD_HPP_)
#define _LSS_HEAT_EULER_SOLVER_METHOD_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heat_euler_coefficients.hpp"
#include "lss_heat_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_utility::container_t;

/**
    explicit_heat_scheme object
 */
class explicit_heat_scheme
{

  public:
    static void rhs(heat_euler_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                    container_t const &input, boundary_1d_pair const &boundary_pair, double const &time,
                    container_t &solution);

    static void rhs_source(heat_euler_coefficients_ptr const &cfs, grid_config_1d_ptr const &grid_config,
                           container_t const &input, container_t const &inhom_input,
                           boundary_1d_pair const &boundary_pair, double const &time, container_t &solution);
};

/**
template <typename fp_type> class heat_euler_solver_method
 object
*/

class heat_euler_solver_method : public heat_explicit_solver_method
{

  private:
    // scheme coefficients:
    heat_euler_coefficients_ptr coefficients_;
    container_t source_;

    explicit heat_euler_solver_method() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heat_euler_solver_method(heat_euler_coefficients_ptr const &coefficients,
                                      grid_config_1d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_euler_solver_method();

    heat_euler_solver_method(heat_euler_solver_method const &) = delete;
    heat_euler_solver_method(heat_euler_solver_method &&) = delete;
    heat_euler_solver_method &operator=(heat_euler_solver_method const &) = delete;
    heat_euler_solver_method &operator=(heat_euler_solver_method &&) = delete;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, container_t &solution) override;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, double const &next_time, std::function<double(double, double)> const &heat_source,
               container_t &solution) override;
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EULER_SOLVER_METHOD_HPP_
