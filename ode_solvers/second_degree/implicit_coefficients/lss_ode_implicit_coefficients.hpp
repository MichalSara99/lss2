/**

    @file      lss_ode_implicit_coefficients.hpp
    @brief     ODE implicit coefficicnts
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_ODE_IMPLICIT_COEFFICIENTS_HPP_)
#define _LSS_ODE_IMPLICIT_COEFFICIENTS_HPP_

#include <functional>

#include "../../../common/lss_range.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../discretization/lss_discretization.hpp"
#include "../../../ode_solvers/lss_ode_discretization_config.hpp"
#include "../../../ode_solvers/transformation/lss_ode_data_transform.hpp"

namespace lss_ode_solvers
{

using lss_utility::range_ptr;
using lss_utility::sptr_t;

/**
    @struct ode_implicit_coefficients
    @brief
**/
struct ode_implicit_coefficients
{
  public:
    // scheme coefficients:
    double lambda_, gamma_;
    std::size_t space_size_;
    range_ptr range_;
    // functional coefficients:
    std::function<double(double)> A_;
    std::function<double(double)> B_;
    std::function<double(double)> C_;

  private:
    void initialize(ode_discretization_config_ptr const &discretization_config);

    void initialize_coefficients(ode_data_transform_ptr const &ode_data_config);

  public:
    ode_implicit_coefficients() = delete;

    explicit ode_implicit_coefficients(ode_data_transform_ptr const &ode_data_config,
                                       ode_discretization_config_ptr const &discretization_config);
};

using ode_implicit_coefficients_ptr = sptr_t<ode_implicit_coefficients>;

} // namespace lss_ode_solvers

#endif ///_LSS_ODE_IMPLICIT_COEFFICIENTS_HPP_
