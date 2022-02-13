/**

    @file      lss_heat_barakat_clark_coefficients.hpp
    @brief     Coefficients for Barakat-Clark heat solvers
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_BARAKAT_CLARK_COEFFICIENTS_HPP_)
#define _LSS_HEAT_BARAKAT_CLARK_COEFFICIENTS_HPP_

#include <functional>

#include "../../../../common/lss_macros.hpp"
#include "../implicit_coefficients/lss_heat_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{
using lss_utility::sptr_t;

struct heat_barakat_clark_coefficients
{
  public:
    // scheme coefficients:
    double k_;
    std::size_t space_size_;
    // functional coefficients:
    std::function<double(double, double)> A_;
    std::function<double(double, double)> B_;
    std::function<double(double, double)> D_;
    std::function<double(double, double)> K_;

  private:
    void initialize_coefficients(heat_coefficients_ptr const &coefficients);

  public:
    heat_barakat_clark_coefficients() = delete;

    explicit heat_barakat_clark_coefficients(heat_coefficients_ptr const &coefficients);
};

using heat_barakat_clark_coefficients_ptr = sptr_t<heat_barakat_clark_coefficients>;

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_BARAKAT_CLARK_COEFFICIENTS_HPP_
