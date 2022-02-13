/**

    @file      lss_heat_coefficients.hpp
    @brief     Coefficients for implicit heat solvers
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_HEAT_COEFFICIENTS_HPP_)
#define _LSS_HEAT_COEFFICIENTS_HPP_

#include <functional>

#include "../../../../common/lss_utility.hpp"
#include "../../../lss_pde_discretization_config.hpp"
#include "../../../transformation/lss_heat_data_transform.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{
using lss_utility::range;
using lss_utility::sptr_t;

struct heat_coefficients
{
  public:
    // scheme coefficients:
    double lambda_, gamma_, delta_, k_;
    std::size_t space_size_;
    range_ptr range_;
    // theta variable:
    double theta_;
    // functional coefficients:
    std::function<double(double, double)> A_;
    std::function<double(double, double)> B_;
    std::function<double(double, double)> D_;

  private:
    void initialize(pde_discretization_config_1d_ptr const &discretization_config);

    void initialize_coefficients(heat_data_transform_1d_ptr const &heat_data_config);

  public:
    heat_coefficients() = delete;

    explicit heat_coefficients(heat_data_transform_1d_ptr const &heat_data_config,
                               pde_discretization_config_1d_ptr const &discretization_config, double const &theta);
};

using heat_coefficients_ptr = sptr_t<heat_coefficients>;

} // namespace one_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HEAT_COEFFICIENTS_HPP_
