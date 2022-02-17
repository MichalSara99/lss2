#include "lss_utility.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>

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

double trapezoidal_method(container_t const &x, container_t const &y)
{
    LSS_ASSERT(x.size() == y.size(), "trapezoidal_method: Sizes of the containers must equal");
    std::size_t n = x.size();
    double sum{0.0};
    for (std::size_t i = 1; i < n; ++i)
    {
        sum += 0.5 * (x[i] - x[i - 1]) * (y[i - 1] + y[i]);
    }
    return sum;
}

heston_exact::heston_exact(double time, double strike, double rate, double lambda, double kappa, double theta,
                           double vol, double rho, double maturity)
    : time_{time}, strike_{strike}, rate_{rate}, lambda_{lambda}, kappa_{kappa}, theta_{theta}, vol_{vol}, rho_{rho},
      maturity_{maturity}
{
    init_integrand();
}

heston_exact::~heston_exact()
{
}

void heston_exact::init_integrand()
{
    integrand_ = [](std::complex<double> K, double X, double vol_0, double tau, double thet, double kappa,
                    double sig_vol, double rho, double gamma) {
        std::complex<double> thetaadj;
        const double omega = kappa * thet;
        const double ksi = sig_vol;
        const double theta = kappa;
        const std::complex<double> t((ksi * ksi) * tau * 0.5, 0.0);
        const std::complex<double> a((2.0 * omega) / (ksi * ksi), 0.0);
        if (gamma == 1.0)
            thetaadj = std::complex<double>(theta, 0.0);
        else
            thetaadj = std::complex<double>(
                (1.0 - gamma) * rho * ksi + sqrt(theta * theta - gamma * (1.0 - gamma) * ksi * ksi), 0.0);

        std::complex<double> im(0.0, 1.0);
        std::complex<double> re(1.0, 0.0);
        std::complex<double> b = (2.0 / (ksi * ksi)) * (thetaadj);
        std::complex<double> c = (K * K - im * K) / (ksi * ksi);
        std::complex<double> d = sqrt(b * b + 4.0 * c);
        std::complex<double> g = (b + d) / 2.0;
        std::complex<double> h = (b + d) / (b - d);
        std::complex<double> f1 = a * (t * g - log((1.0 - h * exp(t * d)) / (1.0 - h)));
        std::complex<double> f2 = g * ((1.0 - exp(t * d)) / (1.0 - h * exp(t * d)));
        std::complex<double> H = exp(f1 + f2 * vol_0);
        // function to be integrated
        std::complex<double> integrand = exp(-im * K * X) * (H / (K * K - im * K));

        return std::real(integrand);
    };
}

void heston_exact::core(double spot, double vol, double time, container_t &x, container_t &y) const
{
    const double ki = 0.5;
    const double tau = maturity_ - time;
    std::complex<double> pass_phi;
    double omega = kappa_ * theta_;
    double ksi = vol_;
    double theta = kappa_;
    const double kmax = ceil(std::max(1000.0, 10.0 / sqrt(vol * tau)));
    const std::size_t n = static_cast<std::size_t>(kmax) * 5;
    container_t int_x(n);
    container_t int_y(n);
    double X = log(spot / strike_) + (rate_)*tau;
    int count = 0;
    for (double phi = 0.000001; phi < kmax; phi += 0.2)
    {
        int_x[count] = phi;
        pass_phi = std::complex<double>(phi, ki);
        int_y[count] = integrand_(pass_phi, X, vol, tau, theta_, kappa_, vol_, rho_, lambda_);
        count += 1;
    }
    x = int_x;
    y = int_y;
}

double heston_exact::call(double spot, double vol) const
{
    container_t x;
    container_t f;
    const double tau = maturity_ - time_;
    core(spot, vol, time_, x, f);
    double price = spot - (1.0 / (pi())) * strike_ * std::exp(-rate_ * tau) * trapezoidal_method(x, f);
    return price;
}

double heston_exact::call(double spot, double vol, double time_to_maturity) const
{
    container_t x;
    container_t f;
    const double tau = time_to_maturity;
    core(spot, vol, time_, x, f);
    double price = spot - (1.0 / (pi())) * strike_ * std::exp(-rate_ * tau) * trapezoidal_method(x, f);
    return price;
}

double heston_exact::put(double spot, double vol) const
{
    container_t x;
    container_t f;
    const double tau = maturity_ - time_;
    core(spot, vol, time_, x, f);
    double price =
        strike_ * exp(-rate_ * tau) - (1.0 / (pi())) * strike_ * exp(-rate_ * tau) * trapezoidal_method(x, f);
    return price;
}

double heston_exact::put(double spot, double vol, double time_to_maturity) const
{
    container_t x;
    container_t f;
    const double tau = time_to_maturity;
    core(spot, vol, time_, x, f);
    double price =
        strike_ * exp(-rate_ * tau) - (1.0 / (pi())) * strike_ * exp(-rate_ * tau) * trapezoidal_method(x, f);
    return price;
}

} // namespace lss_utility
