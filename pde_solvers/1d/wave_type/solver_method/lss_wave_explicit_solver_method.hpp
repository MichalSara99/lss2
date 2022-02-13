/**

    @file      lss_wave_explicit_solver_method.hpp
    @brief     Abstract explicit solver method for wave problems
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_WAVE_EXPLICIT_SOLVER_METHOD_HPP_)
#define _LSS_WAVE_EXPLICIT_SOLVER_METHOD_HPP_

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
wave_euler_solver_method object
*/
class wave_explicit_solver_method
{
  protected:
    grid_config_1d_ptr grid_cfg_;

    explicit wave_explicit_solver_method() = delete;

  public:
    explicit wave_explicit_solver_method(grid_config_1d_ptr const &grid_config);

    virtual ~wave_explicit_solver_method();

    virtual void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                               boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                               container_t &solution) = 0;

    virtual void solve_initial(container_t &prev_solution_0, container_t &prev_solution_1,
                               boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                               std::function<double(double, double)> const &wave_source, container_t &solution) = 0;

    virtual void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                                container_t &solution) = 0;

    virtual void solve_terminal(container_t &prev_solution_0, container_t &prev_solution_1,
                                boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                                std::function<double(double, double)> const &wave_source, container_t &solution) = 0;

    virtual void solve(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       container_t &solution) = 0;

    virtual void solve(container_t &prev_solution_0, container_t &prev_solution_1,
                       boundary_1d_pair const &boundary_pair, double const &time, double const &next_time,
                       std::function<double(double, double)> const &wave_source, container_t &solution) = 0;
};

using wave_explicit_solver_method_ptr = sptr_t<wave_explicit_solver_method>;

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_WAVE_EXPLICIT_SOLVER_METHOD_HPP_
