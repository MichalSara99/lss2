#include "lss_wave_implicit_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

void wave_implicit_coefficients::initialize(pde_discretization_config_1d_ptr const &discretization_config)
{
    // get space range:
    range_ = discretization_config->space_range();
    // get time step:
    k_ = discretization_config->time_step();
    // size of spaces discretization:
    space_size_ = discretization_config->number_of_space_points();
    const double one = 1.0;
    const double two = 2.0;
    const double h = one / (space_size_ - 1);
    // calculate scheme coefficients:
    lambda_ = one / (k_ * k_);
    gamma_ = one / (two * k_);
    delta_ = one / (h * h);
    rho_ = one / (two * h);
}

void wave_implicit_coefficients::initialize_coefficients(wave_data_transform_1d_ptr const &wave_data_config)
{
    // save coefficients locally:
    auto const a = wave_data_config->a_coefficient();
    auto const b = wave_data_config->b_coefficient();
    auto const c = wave_data_config->c_coefficient();
    auto const d = wave_data_config->d_coefficient();

    const double half = 0.5;
    const double quater = 0.25;

    A_ = [=](double t, double x) { return quater * (delta_ * b(t, x) - rho_ * c(t, x)); };
    B_ = [=](double t, double x) { return quater * (delta_ * b(t, x) + rho_ * c(t, x)); };
    C_ = [=](double t, double x) { return half * (delta_ * b(t, x) - half * d(t, x)); };
    D_ = [=](double t, double x) { return (lambda_ - gamma_ * a(t, x)); };
    E_ = [=](double t, double x) { return (lambda_ + gamma_ * a(t, x)); };
}

wave_implicit_coefficients::wave_implicit_coefficients(wave_data_transform_1d_ptr const &wave_data_config,
                                                       pde_discretization_config_1d_ptr const &discretization_config)
{
    initialize(discretization_config);
    initialize_coefficients(wave_data_config);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
