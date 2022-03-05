/**

    @file      lss_heat_splitting_method_3d.hpp
    @brief     Abstract class for splitting methods for 3D heat problems
    @details   ~
    @author    Michal Sara
    @date      5.03.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_SPLITTING_METHOD_3D_HPP_)
#define _LSS_HEAT_SPLITTING_METHOD_3D_HPP_

#include <functional>

#include "../../../../../../boundaries/lss_boundary.hpp"
#include "../../../../../../common/lss_enumerations.hpp"
#include "../../../../../../common/lss_macros.hpp"
#include "../../../../../../common/lss_utility.hpp"
#include "../../../../../../containers/lss_matrix_3d.hpp"

namespace lss_pde_solvers
{

namespace three_dimensional
{
using lss_boundary::boundary_3d_pair;
using lss_boundary::boundary_3d_ptr;
using lss_containers::matrix_3d;
using lss_utility::sptr_t;

/**
    heat_splitting_method_3d object
 */
class heat_splitting_method_3d
{
  public:
    explicit heat_splitting_method_3d();

    ~heat_splitting_method_3d();

    virtual void solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
                       boundary_3d_pair const &y_boundary_pair, boundary_3d_pair const &z_boundary_pair,
                       double const &time, matrix_3d &solution) = 0;

    virtual void solve(matrix_3d const &prev_solution, boundary_3d_pair const &x_boundary_pair,
                       boundary_3d_pair const &y_boundary_pair, boundary_3d_pair const &z_boundary_pair,
                       double const &time, std::function<double(double, double, double)> const &heat_source,
                       matrix_3d &solution) = 0;
};

using heat_splitting_method_3d_ptr = sptr_t<heat_splitting_method_3d>;

} // namespace three_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_SPLITTING_METHOD_3D_HPP_
