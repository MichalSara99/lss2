#include "lss_utility.hpp"

#include <cmath>
#include <limits>

namespace lss_utility
{

constexpr double NaN()
{
    return std::numeric_limits<double>::quiet_NaN();
}

void valcopy(container_t &dst, container_t const &src, std::size_t dst_offset)
{
    LSS_ASSERT((dst.size() - src.size()) >= 0, "src must be at least al long as dst");
    LSS_ASSERT(dst_offset >= 0, "dst_offset must be non-negative");
    for (std::size_t t = 0; t < src.size(); ++t)
    {
        dst[t + dst_offset] = src[t];
    }
}

void copy(container_t &dst, thrust::host_vector<double> const &src)
{
    for (std::size_t t = 0; t < dst.size(); ++t)
    {
        dst[t] = src[t];
    }
}

void copy(thrust::host_vector<double> &dst, container_t const &src)
{
    for (std::size_t t = 0; t < dst.size(); ++t)
    {
        dst[t] = src[t];
    }
}

double norm_cdf(double x)
{
    double ind = double{};
    if (x <= 0.0)
        ind = 1.0;
    x = std::abs(x);
    double const cst = 1.0 / (std::sqrt(2.0 * pi()));
    double const first = std::exp(-0.5 * x * x);
    double const second = 0.226 + 0.64 * x + 0.33 * std::sqrt(x * x + 3.0);
    double const res = 1.0 - ((first / second) * cst);
    return std::abs(ind - res);
}

black_scholes_exact::black_scholes_exact()
{
}

black_scholes_exact::black_scholes_exact(double time, double strike, double rate, double volatility, double maturity)
    : time_{time}, strike_{strike}, rate_{rate}, vol_{volatility}, maturity_{maturity}
{
}

black_scholes_exact ::~black_scholes_exact()
{
}

double black_scholes_exact::call(double spot) const
{
    double const tau = maturity_ - time_;
    double const s_tau = std::sqrt(tau);
    double const d_1 = (std::log(spot / strike_) + (rate_ + 0.5 * vol_ * vol_) * tau) / (vol_ * s_tau);
    double const d_2 = d_1 - vol_ * s_tau;
    double const result = norm_cdf(d_1) * spot - (norm_cdf(d_2) * strike_ * std::exp(-rate_ * tau));
    return result;
}

double black_scholes_exact::call(double spot, double time_to_maturity) const
{
    double const tau = time_to_maturity;
    double const s_tau = std::sqrt(tau);
    double const d_1 = (std::log(spot / strike_) + (rate_ + 0.5 * vol_ * vol_) * tau) / (vol_ * s_tau);
    double const d_2 = d_1 - vol_ * s_tau;
    double const result = norm_cdf(d_1) * spot - (norm_cdf(d_2) * strike_ * std::exp(-rate_ * tau));
    return result;
}

double black_scholes_exact::put(double spot) const
{
    double const call_p = call(spot);
    double const tau = maturity_ - time_;
    return (strike_ * std::exp(-rate_ * tau) - spot + call_p);
}

double black_scholes_exact::put(double spot, double time_to_maturity) const
{
    double const call_p = call(spot, time_to_maturity);
    double const tau = time_to_maturity;
    return (strike_ * std::exp(-rate_ * tau) - spot + call_p);
}

} // namespace lss_utility
