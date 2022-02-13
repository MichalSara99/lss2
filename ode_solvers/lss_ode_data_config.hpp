/**

    @file      lss_ode_data_config.hpp
    @brief     Represents data configuartion for ODE solvers
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_ODE_DATA_CONFIG_HPP_)
#define _LSS_ODE_DATA_CONFIG_HPP_

#include <functional>

#include "../common/lss_macros.hpp"
#include "../common/lss_utility.hpp"

namespace lss_ode_solvers
{

using lss_utility::sptr_t;

/**
    @struct ode_coefficient_data_config
    @brief
**/
struct ode_coefficient_data_config
{
  private:
    std::function<double(double)> a_coeff_;
    std::function<double(double)> b_coeff_;

    explicit ode_coefficient_data_config() = delete;

    void initialize();

  public:
    explicit ode_coefficient_data_config(std::function<double(double)> const &a_coefficient,
                                         std::function<double(double)> const &b_coefficient);

    /**
        @brief First-derivative function coefficient
        @retval function
    **/
    LSS_API std::function<double(double)> const &a_coefficient() const;

    /**
        @brief function coefficient
        @retval function
    **/
    LSS_API std::function<double(double)> const &b_coefficient() const;
};

using ode_coefficient_data_config_ptr = sptr_t<ode_coefficient_data_config>;

/**
    @struct ode_nonhom_data_config
    @brief
**/
struct ode_nonhom_data_config
{
  private:
    std::function<double(double)> nonhom_fun_;

    explicit ode_nonhom_data_config() = delete;

  public:
    explicit ode_nonhom_data_config(std::function<double(double)> const &nonhom_fun);

    /**
        @brief nonhomogeneous function
        @retval function
    **/
    LSS_API std::function<double(double)> const &nonhom_function() const;
};

using ode_nonhom_data_config_ptr = sptr_t<ode_nonhom_data_config>;

/**
    @struct ode_data_config
    @brief
**/
struct ode_data_config
{
  private:
    ode_coefficient_data_config_ptr coefficient_data_cfg_;
    ode_nonhom_data_config_ptr nonhom_data_cfg_;

    void initialize();

    explicit ode_data_config() = delete;

  public:
    explicit ode_data_config(ode_coefficient_data_config_ptr const &coefficient_data_config,
                             ode_nonhom_data_config_ptr const &nonhom_data_config = nullptr);

    ~ode_data_config();

    /**
        @brief Configuration object for nonhomogeneous function
        @retval ode_nonhom_data_config_ptr
    **/
    LSS_API ode_nonhom_data_config_ptr const &nonhom_data_config() const;

    /**
        @brief First-derivative function coefficient
        @retval function
    **/
    LSS_API std::function<double(double)> const &a_coefficient() const;

    /**
        @brief function coefficient
        @retval function
    **/
    LSS_API std::function<double(double)> const &b_coefficient() const;
};

using ode_data_config_ptr = sptr_t<ode_data_config>;

} // namespace lss_ode_solvers

#endif ///_LSS_ODE_DATA_CONFIG_HPP_
