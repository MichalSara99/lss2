/**

    @file      lss_heat_saulyev_coefficients.hpp
    @brief     Coefficients for Saulyev heat solvers
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_SAULYEV_COEFFICIENTS_HPP_)
#define _LSS_HEAT_SAULYEV_COEFFICIENTS_HPP_

#include <functional>

#include "../../../../common/lss_macros.hpp"
#include "../implicit_coefficients/lss_heat_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{
using lss_utility::sptr_t;

struct heat_saulyev_coefficients
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
    heat_saulyev_coefficients() = delete;

    explicit heat_saulyev_coefficients(heat_coefficients_ptr const &coefficients);
};

using heat_saulyev_coefficients_ptr = sptr_t<heat_saulyev_coefficients>;

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_SAULYEV_COEFFICIENTS_HPP_
