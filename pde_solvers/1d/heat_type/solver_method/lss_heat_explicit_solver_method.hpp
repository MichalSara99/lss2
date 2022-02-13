/**

    @file      lss_heat_explicit_solver_method.hpp
    @brief     Abstract Euler solver method for heat problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EXPLICIT_SOLVER_METHOD_HPP_)
#define _LSS_HEAT_EXPLICIT_SOLVER_METHOD_HPP_

#include "../../../../boundaries/lss_boundary.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../discretization/lss_grid_config.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

using lss_boundary::boundary_1d_pair;
using lss_grids::grid_config_1d_ptr;
using lss_utility::container_t;
using lss_utility::sptr_t;

/**
 heat_explicit_solver_method  object

*/
class heat_explicit_solver_method
{

  protected:
    grid_config_1d_ptr grid_cfg_;

    explicit heat_explicit_solver_method() = delete;

  public:
    explicit heat_explicit_solver_method(grid_config_1d_ptr const &grid_config);

    virtual ~heat_explicit_solver_method();

    virtual void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
                       double const &time, container_t &solution) = 0;

    virtual void solve(container_t &prev_solution, boundary_1d_pair const &boundary_pair, std::size_t const &time_idx,
                       double const &time, double const &next_time,
                       std::function<double(double, double)> const &heat_source, container_t &solution) = 0;
};

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EXPLICIT_SOLVER_METHOD_HPP_
