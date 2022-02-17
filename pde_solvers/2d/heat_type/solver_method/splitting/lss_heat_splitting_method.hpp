/**

    @file      lss_heat_splitting_method.hpp
    @brief     Abstract class for splitting methods for 2D heat problems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_SPLITTING_METHOD_HPP_)
#define _LSS_HEAT_SPLITTING_METHOD_HPP_

#include <functional>

#include "../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../common/lss_enumerations.hpp"
#include "../../../../../common/lss_macros.hpp"
#include "../../../../../common/lss_utility.hpp"
#include "../../../../../containers/lss_matrix_2d.hpp"

namespace lss_pde_solvers
{

namespace two_dimensional
{
using lss_boundary::boundary_2d_pair;
using lss_containers::matrix_2d;
using lss_utility::sptr_t;

/**
    heat_splitting_method object
 */
class heat_splitting_method
{
  public:
    explicit heat_splitting_method();

    ~heat_splitting_method();

    virtual void solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizontal_boundary_pair,
                       boundary_2d_pair const &vertical_boundary_pair, double const &time, matrix_2d &solution) = 0;

    virtual void solve(matrix_2d const &prev_solution, boundary_2d_pair const &horizontal_boundary_pair,
                       boundary_2d_pair const &vertical_boundary_pair, double const &time,
                       std::function<double(double, double)> const &heat_source, matrix_2d &solution) = 0;
};

using heat_splitting_method_ptr = sptr_t<heat_splitting_method>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_SPLITTING_METHOD_HPP_
