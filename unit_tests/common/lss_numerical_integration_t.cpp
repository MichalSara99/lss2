#include "lss_numerical_integration_t.hpp"

#include <functional>

using lss_utility::container_t;
using lss_utility::trapezoidal_method;

void trapezoid_integral_1_test()
{
    auto const integrand = [](double x) { return x * x; };
    auto const low = 0.0;
    auto const high = 1.0;
    const std::size_t n = 1000;
    const double step = 1.0 / static_cast<double>(n - 1);

    container_t x(low, n);
    container_t f(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        x[i] = static_cast<double>(low + i * step);
        f[i] = integrand(x[i]);
    }
    auto const true_value = 1.0 / 3.;
    auto const result = trapezoidal_method(x, f);
    LSS_ASSERT(std::abs(true_value - result) < 1e-6, "Integral value must be approximately same");
}

void trapezoid_integral_2_test()
{
}
