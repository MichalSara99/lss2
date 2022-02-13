/**

    @file      lss_utility.hpp
    @brief     Common utilities
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_UTILITY_HPP_)
#define _LSS_UTILITY_HPP_

#include <memory>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <valarray>

#include "lss_macros.hpp"

namespace lss_utility
{

using container_t = std::valarray<double>;

void valcopy(container_t &dst, container_t const &src, std::size_t dst_offset = 0);

void copy(container_t &dst, thrust::host_vector<double> const &src);

void copy(thrust::host_vector<double> &dst, container_t const &src);

/**
    @brief  Pi constant
    @retval
**/
static constexpr double pi()
{
    return 3.14159265358979;
}

/**
    @brief  Quiet NaN
    @retval
**/
static constexpr double NaN();

template <typename T> using sptr_t = std::shared_ptr<T>;

template <typename T> using uptr_t = std::unique_ptr<T>;

/**
    @brief  Cumulative distribution function for normal random variate
    @param  x - value
    @retval
**/
double norm_cdf(double x);

/**
    @struct black_scholes_exact
    @brief  represents exact BS formula
**/
struct black_scholes_exact
{
  private:
    double time_;
    double strike_;
    double vol_;
    double maturity_;
    double rate_;

  protected:
    explicit black_scholes_exact();

  public:
    /**
        @brief black_scholes_exact object constructor
        @param time - value time
        @param strike - strie value
        @param rate - rate value
        @param volatility - volatility value
        @param maturity - maturity value
    **/
    explicit black_scholes_exact(double time, double strike, double rate, double volatility, double maturity);

    ~black_scholes_exact();

    /**
        @brief  value of call
        @param  spot - spot value
        @retval
    **/
    double call(double spot) const;

    /**
        @brief  value of call
        @param  spot - spot value
        @param  time_to_maturity - time to maturity
        @retval
    **/
    double call(double spot, double time_to_maturity) const;

    /**
        @brief  value of put
        @param  spot - spot value
        @retval
    **/
    double put(double spot) const;

    /**
        @brief  value of put
        @param  spot - spot value
        @param  time_to_maturity - time to maturity
        @retval
    **/
    double put(double spot, double time_to_maturity) const;
};

} // namespace lss_utility

#endif ///_LSS_UTILITY_HPP_
