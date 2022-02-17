#include "lss_heat_coefficients_2d.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

void heat_coefficients_2d::initialize(pde_discretization_config_2d_ptr const &discretization_config,
                                      splitting_method_config_ptr const &splitting_config)
{
    // get space ranges:
    const auto &spaces = discretization_config->space_range();
    // across X:
    rangex_ = std::get<0>(spaces);
    // across Y:
    rangey_ = std::get<1>(spaces);
    // time step:
    const double k = discretization_config->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_config->number_of_space_points();
    space_size_x_ = std::get<0>(space_sizes);
    space_size_y_ = std::get<1>(space_sizes);
    // calculate scheme coefficients:
    const double one = 1.0;
    const double half = 0.5;
    const double quarter = 0.25;
    auto const &h_1 = one / (space_size_x_ - 1);
    auto const &h_2 = one / (space_size_y_ - 1);

    k_ = k;
    alpha_ = k_ / (h_1 * h_1);
    beta_ = k_ / (h_2 * h_2);
    gamma_ = quarter * k_ / (h_1 * h_2);
    delta_ = half * k_ / h_1;
    ni_ = half * k_ / h_2;
    rho_ = k;
    if (splitting_config != nullptr)
        zeta_ = splitting_config->weighting_value();
}

void heat_coefficients_2d::initialize_coefficients(heat_data_transform_2d_ptr const &heat_data_config)
{
    // save coefficients locally:
    auto const a = heat_data_config->a_coefficient();
    auto const b = heat_data_config->b_coefficient();
    auto const c = heat_data_config->c_coefficient();
    auto const d = heat_data_config->d_coefficient();
    auto const e = heat_data_config->e_coefficient();
    auto const f = heat_data_config->f_coefficient();

    const double two = 2.0;
    const double half = 0.5;

    M_ = [=](double t, double x, double y) { return (alpha_ * a(t, x, y) - delta_ * d(t, x, y)); };
    M_tilde_ = [=](double t, double x, double y) { return (beta_ * b(t, x, y) - ni_ * e(t, x, y)); };
    P_ = [=](double t, double x, double y) { return (alpha_ * a(t, x, y) + delta_ * d(t, x, y)); };
    P_tilde_ = [=](double t, double x, double y) { return (beta_ * b(t, x, y) + ni_ * e(t, x, y)); };
    Z_ = [=](double t, double x, double y) { return (two * alpha_ * a(t, x, y) - half * rho_ * f(t, x, y)); };
    W_ = [=](double t, double x, double y) { return (two * beta_ * b(t, x, y) - half * rho_ * f(t, x, y)); };
    C_ = [=](double t, double x, double y) { return c(t, x, y); };
    D_ = [=](double t, double x, double y) { return d(t, x, y); };
    E_ = [=](double t, double x, double y) { return e(t, x, y); };
    F_ = [=](double t, double x, double y) { return f(t, x, y); };
}

heat_coefficients_2d::heat_coefficients_2d(heat_data_transform_2d_ptr const &heat_data_config,
                                           pde_discretization_config_2d_ptr const &discretization_config,
                                           splitting_method_config_ptr const splitting_config, double const &theta)
    : theta_{theta}
{
    initialize(discretization_config, splitting_config);
    initialize_coefficients(heat_data_config);
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
