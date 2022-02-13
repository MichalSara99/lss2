/**

    @file      lss_ode_discretization_config.hpp
    @brief     ODE discratization configuration
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_ODE_DISCRETIZATION_CONFIG_HPP_)
#define _LSS_ODE_DISCRETIZATION_CONFIG_HPP_

#include "../common/lss_range.hpp"
#include "../common/lss_utility.hpp"
#include "../discretization/lss_discretization_config.hpp"

namespace lss_ode_solvers
{

using lss_utility::range_ptr;
using lss_utility::sptr_t;

/**
    ode_discretization_config structure
 */
struct ode_discretization_config final : public lss_discretization::discretization_config_1d
{
  private:
    explicit ode_discretization_config() = delete;

  public:
    explicit ode_discretization_config(range_ptr const &space_range, std::size_t const &number_of_space_points);

    ~ode_discretization_config();
};

using ode_discretization_config_ptr = sptr_t<ode_discretization_config>;

} // namespace lss_ode_solvers

#endif ///_LSS_ODE_DISCRETIZATION_CONFIG_HPP_
