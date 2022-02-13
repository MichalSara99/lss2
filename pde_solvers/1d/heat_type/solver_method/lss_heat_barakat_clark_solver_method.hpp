/**

    @file      lss_heat_barakat_clark_solver_method.hpp
    @brief     Barakat-Clark solver method for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_BARAKAT_CLARK_SOLVER_METHOD_HPP_)
#define _LSS_HEAT_BARAKAT_CLARK_SOLVER_METHOD_HPP_

#include <future>

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heat_barakat_clark_coefficients.hpp"
#include "lss_heat_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_utility::container_t;
using lss_utility::sptr_t;

/**
 heat_barakat_clark_solver_method  object

*/

class heat_barakat_clark_solver_method : public heat_explicit_solver_method
{
    typedef std::function<void(container_t &, container_t const &, double, double)> sweeper_fun;

  private:
    // constant coeffs:
    const double cone_ = 1.0;
    const double chalf_ = 0.5;
    const double czero_ = 0.0;
    // scheme coefficients:
    heat_barakat_clark_coefficients_ptr coefficients_;
    // sweepers:
    sweeper_fun up_sweeper_, down_sweeper_;
    // containers:
    container_t source_dummy_, source_, source_next_;

    explicit heat_barakat_clark_solver_method() = delete;

    void initialize(bool is_heat_sourse_set);

  public:
    explicit heat_barakat_clark_solver_method(heat_barakat_clark_coefficients_ptr const &coefficients,
                                              grid_config_1d_ptr const &grid_config, bool is_heat_sourse_set);

    ~heat_barakat_clark_solver_method();

    heat_barakat_clark_solver_method(heat_barakat_clark_solver_method const &) = delete;
    heat_barakat_clark_solver_method(heat_barakat_clark_solver_method &&) = delete;
    heat_barakat_clark_solver_method &operator=(heat_barakat_clark_solver_method const &) = delete;
    heat_barakat_clark_solver_method &operator=(heat_barakat_clark_solver_method &&) = delete;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, container_t &solution) override;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, double const &next_time, std::function<double(double, double)> const &heat_source,
               container_t &solution) override;
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_BARAKAT_CLARK_SOLVER_METHOD_HPP_
