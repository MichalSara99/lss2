/**

    @file      lss_heat_saulyev_solver_method.hpp
    @brief     Saulyev solver method for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_SAULYEV_SOLVER_METHOD_HPP_)
#define _LSS_HEAT_SAULYEV_SOLVER_METHOD_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../discretization/lss_discretization.hpp"
#include "../../../../discretization/lss_grid_config.hpp"
#include "../explicit_coefficients/lss_heat_saulyev_coefficients.hpp"
#include "lss_heat_explicit_solver_method.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_1d_ptr;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;
using lss_utility::sptr_t;

/**
 heat_saulyev_solver_method  object

*/

class heat_saulyev_solver_method : public heat_explicit_solver_method
{

    typedef std::function<void(container_t &, container_t const &, double, double)> sweeper_fun;

  private:
    // constant coeffs:
    const double czero_ = 0.0;
    const double cone_ = 1.0;
    // scheme coefficients:
    heat_saulyev_coefficients_ptr coefficients_;
    // sweepers:
    sweeper_fun up_sweeper_, down_sweeper_;
    // containers for source:
    container_t source_dummy_, source_, source_next_;

    explicit heat_saulyev_solver_method() = delete;

    void initialize(bool is_heat_source_set);

  public:
    explicit heat_saulyev_solver_method(heat_saulyev_coefficients_ptr const &coefficients,
                                        grid_config_1d_ptr const &grid_config, bool is_heat_source_set);

    ~heat_saulyev_solver_method();

    heat_saulyev_solver_method(heat_saulyev_solver_method const &) = delete;
    heat_saulyev_solver_method(heat_saulyev_solver_method &&) = delete;
    heat_saulyev_solver_method &operator=(heat_saulyev_solver_method const &) = delete;
    heat_saulyev_solver_method &operator=(heat_saulyev_solver_method &&) = delete;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, container_t &solution) override;

    void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
               double const &time, double const &next_time, std::function<double(double, double)> const &heat_source,
               container_t &solution) override;
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_SAULYEV_SOLVER_METHOD_HPP_
