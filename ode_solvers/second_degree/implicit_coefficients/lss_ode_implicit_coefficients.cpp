#include "lss_ode_implicit_coefficients.hpp"

namespace lss_ode_solvers
{

ode_implicit_coefficients::ode_implicit_coefficients(ode_data_transform_ptr const &ode_data_config,
                                                     ode_discretization_config_ptr const &discretization_config)
{
    initialize(discretization_config);
    initialize_coefficients(ode_data_config);
}

void ode_implicit_coefficients::initialize_coefficients(ode_data_transform_ptr const &ode_data_config)
{
    // save coefficients locally:
    auto const a = ode_data_config->a_coefficient();
    auto const b = ode_data_config->b_coefficient();
    const double two = static_cast<double>(2.0);
    A_ = [=](double x) { return (lambda_ - gamma_ * a(x)); };
    B_ = [=](double x) { return (lambda_ + gamma_ * a(x)); };
    C_ = [=](double x) { return (b(x) - two * lambda_); };
}

void ode_implicit_coefficients::initialize(ode_discretization_config_ptr const &discretization_config)
{
    // get space range:
    range_ = discretization_config->space_range();
    // size of spaces discretization:
    space_size_ = discretization_config->number_of_space_points();
    const double one = static_cast<double>(1.0);
    const double two = static_cast<double>(2.0);
    const double h = one / (space_size_ - 1);
    // calculate scheme coefficients:
    lambda_ = one / (h * h);
    gamma_ = one / (two * h);
}

} // namespace lss_ode_solvers
