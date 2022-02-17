/**

    @file      lss_heat_explicit_solver_method_2d.hpp
    @brief     Abstract explicit Euler method for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EXPLICIT_SOLVER_METHOD_2D_HPP_)
#define _LSS_HEAT_EXPLICIT_SOLVER_METHOD_2D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../containers/lss_matrix_2d.hpp"
#include "../../../../discretization/lss_grid_config.hpp"

namespace lss_pde_solvers
{

using lss_containers::matrix_2d;
using lss_grids::grid_config_2d_ptr;
using lss_utility::sptr_t;

namespace two_dimensional
{

class heat_explicit_solver_method_2d
{

  protected:
    grid_config_2d_ptr grid_cfg_;

    explicit heat_explicit_solver_method_2d() = delete;

  public:
    explicit heat_explicit_solver_method_2d(grid_config_2d_ptr const &grid_config);

    virtual ~heat_explicit_solver_method_2d();

    virtual void solve(matrix_2d &prev_solution, double const &time, matrix_2d &solution) = 0;

    virtual void solve(matrix_2d &prev_solution, double const &time,
                       std::function<double(double, double, double)> const &heat_source, matrix_2d &solution) = 0;
};

using heat_explicit_solver_method_2d_ptr = sptr_t<heat_explicit_solver_method_2d>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EXPLICIT_SOLVER_METHOD_2D_HPP_
