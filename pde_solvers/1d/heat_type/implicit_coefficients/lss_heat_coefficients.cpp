#include "lss_heat_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

void heat_coefficients::initialize(pde_discretization_config_1d_ptr const &discretization_config)
{
    // get space range:
    range_ = discretization_config->space_range();
    // get time step:
    k_ = discretization_config->time_step();
    // size of spaces discretization:
    space_size_ = discretization_config->number_of_space_points();
    const double one = 1.0;
    const double two = 2.0;
    const double half = 0.5;
    auto const h = one / (space_size_ - 1);
    // calculate scheme coefficients:
    lambda_ = k_ / (h * h);
    gamma_ = k_ / (two * h);
    delta_ = half * k_;
}

void heat_coefficients::initialize_coefficients(heat_data_transform_1d_ptr const &heat_data_config)
{
    // save coefficients locally:
    auto const a = heat_data_config->a_coefficient();
    auto const b = heat_data_config->b_coefficient();
    auto const c = heat_data_config->c_coefficient();

    A_ = [=](double t, double x) { return (lambda_ * a(t, x) - gamma_ * b(t, x)); };
    B_ = [=](double t, double x) { return (lambda_ * a(t, x) - delta_ * c(t, x)); };
    D_ = [=](double t, double x) { return (lambda_ * a(t, x) + gamma_ * b(t, x)); };
}

heat_coefficients::heat_coefficients(heat_data_transform_1d_ptr const &heat_data_config,
                                     pde_discretization_config_1d_ptr const &discretization_config, double const &theta)
    : theta_{theta}
{
    initialize(discretization_config);
    initialize_coefficients(heat_data_config);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
