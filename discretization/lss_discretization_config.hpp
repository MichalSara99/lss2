/**

    @file      lss_discretization_config.hpp
    @brief     Represents general 1D discretization configuration
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once

#if !defined(_LSS_DISCRETIZATION_CONFIG_HPP_)
#define _LSS_DISCRETIZATION_CONFIG_HPP_

#include "../common/lss_macros.hpp"
#include "../common/lss_range.hpp"
#include "../common/lss_utility.hpp"

namespace lss_discretization
{

using lss_utility::range_ptr;
using lss_utility::sptr_t;

/**
    @struct discretization_config_1d
    @brief discretization_config structure
**/
struct discretization_config_1d
{
  protected:
    range_ptr space_range_;
    std::size_t number_of_space_points_;

    explicit discretization_config_1d() = delete;

  public:
    explicit discretization_config_1d(range_ptr const &space_range, std::size_t const &number_of_space_points);

    ~discretization_config_1d();

    /**
        @brief Space range
        @retval
    **/
    LSS_API range_ptr const &space_range() const;

    /**
        @brief Number of space points in discretization
        @retval
    **/
    LSS_API std::size_t number_of_space_points() const;

    /**
        @brief  Value for space step
        @retval
    **/
    LSS_API double space_step() const;
};

using discretization_config_1d_ptr = sptr_t<discretization_config_1d>;

} // namespace lss_discretization

#endif ///_LSS_DISCRETIZATION_CONFIG_HPP_
