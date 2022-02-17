/**

    @file      lss_heat_euler_coefficients_2d.hpp
    @brief     Explicit Euler coefficients for 2D heat porblems
    @details   ~
    @author    Michal Sara
    @date      14.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_EULER_COEFFICIENTS_2D_HPP_)
#define _LSS_HEAT_EULER_COEFFICIENTS_2D_HPP_

#include "../../../../common/lss_utility.hpp"
#include "../implicit_coefficients/lss_heat_coefficients_2d.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_utility::sptr_t;

/**
    heat_euler_coefficients_2d object
 */
struct heat_euler_coefficients_2d
{
  public:
    // scheme constant coefficients:
    double rho_, k_;
    std::size_t space_size_x_, space_size_y_;
    // functional coefficients:
    std::function<double(double, double, double)> M_;
    std::function<double(double, double, double)> M_tilde_;
    std::function<double(double, double, double)> P_;
    std::function<double(double, double, double)> P_tilde_;
    std::function<double(double, double, double)> Z_;
    std::function<double(double, double, double)> W_;
    std::function<double(double, double, double)> C_;

  private:
    void initialize_coefficients(heat_coefficients_2d_ptr const &coefficients);

  public:
    heat_euler_coefficients_2d() = delete;

    heat_euler_coefficients_2d(heat_coefficients_2d_ptr const &coefficients);
};

using heat_euler_coefficients_2d_ptr = sptr_t<heat_euler_coefficients_2d>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_EULER_COEFFICIENTS_2D_HPP_
