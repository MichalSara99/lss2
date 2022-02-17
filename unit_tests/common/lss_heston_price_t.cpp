#include "lss_heston_price_t.hpp"

#include <functional>

using lss_utility::container_t;
using lss_utility::heston_exact;

void heston_price_test()
{
    auto const spot = 100.0;
    auto const strike = 100.0;
    auto const rate = 0.05;
    auto const t = 0.0;
    auto const mat = 0.5;
    auto const rho = 0.0;
    auto const kappa = 2.0;
    auto const theta = 0.01;
    auto const lambda = 0.0;
    auto const vol_vol = 0.225;
    auto const vol = 0.01;

    heston_exact heston(t, strike, rate, lambda, kappa, theta, vol_vol, rho, mat);
    auto const call = heston.call(spot, vol);
    auto const put = heston.put(spot, vol);
    auto const true_call = 4.0852;
    auto const true_put = 1.6162;

    LSS_ASSERT(std::abs(true_call - call) < 1e-4, "Integral value must be approximately same");
    LSS_ASSERT(std::abs(true_put - put) < 1e-4, "Integral value must be approximately same");
}
