#include "lss_hhw_implicit_coefficients.hpp"

namespace lss_pde_solvers
{
namespace three_dimensional
{

void hhw_implicit_coefficients::initialize(pde_discretization_config_3d_ptr const &discretization_config,
                                           splitting_method_config_ptr const &splitting_config)
{
    // get space ranges:
    const auto &spaces = discretization_config->space_range();
    // across X:
    rangex_ = std::get<0>(spaces);
    // across Y:
    rangey_ = std::get<1>(spaces);
    // across Z:
    rangez_ = std::get<2>(spaces);
    // time step:
    const double k = discretization_config->time_step();
    // size of spaces discretization:
    const auto &space_sizes = discretization_config->number_of_space_points();
    space_size_x_ = std::get<0>(space_sizes);
    space_size_y_ = std::get<1>(space_sizes);
    space_size_z_ = std::get<2>(space_sizes);
    // calculate scheme coefficients:
    const double one = 1.0;
    const double half = 0.5;
    const double quarter = 0.25;
    auto const h_1 = one / (space_size_x_ - 1);
    auto const h_2 = one / (space_size_y_ - 1);
    auto const h_3 = one / (space_size_z_ - 1);

    k_ = k;
    alpha_1_ = k_ / (h_1 * h_1);
    alpha_2_ = k_ / (h_2 * h_2);
    alpha_3_ = k_ / (h_3 * h_3);
    beta_1_ = quarter * k_ / (h_1 * h_2);
    beta_2_ = quarter * k_ / (h_1 * h_3);
    beta_3_ = quarter * k_ / (h_2 * h_3);
    gamma_1_ = half * k_ / (h_1);
    gamma_2_ = half * k_ / (h_2);
    gamma_3_ = half * k_ / (h_3);
    rho_ = k;
    if (splitting_config != nullptr)
        zeta_ = splitting_config->weighting_value();
}

void hhw_implicit_coefficients::initialize_coefficients(heat_data_transform_3d_ptr const &heat_data_config)
{
    // save coefficients locally:
    auto const a = heat_data_config->a_coefficient();
    auto const b = heat_data_config->b_coefficient();
    auto const c = heat_data_config->c_coefficient();
    auto const d = heat_data_config->d_coefficient();
    auto const e = heat_data_config->e_coefficient();
    auto const f = heat_data_config->f_coefficient();
    auto const g = heat_data_config->g_coefficient();
    auto const h = heat_data_config->h_coefficient();
    auto const i = heat_data_config->i_coefficient();
    auto const j = heat_data_config->j_coefficient();

    const double two = 2.0;
    const double half = 0.5;
    const double sixth = 1.0 / 6.0;

    M_1_ = [=](double t, double x, double y, double z) {
        return (alpha_1_ * a(t, x, y, z) - gamma_1_ * g(t, x, y, z));
    };
    M_2_ = [=](double t, double x, double y, double z) {
        return (alpha_2_ * b(t, x, y, z) - gamma_2_ * h(t, x, y, z));
    };
    M_3_ = [=](double t, double x, double y, double z) {
        return (alpha_3_ * c(t, x, y, z) - gamma_3_ * i(t, x, y, z));
    };
    P_1_ = [=](double t, double x, double y, double z) {
        return (alpha_1_ * a(t, x, y, z) + gamma_1_ * g(t, x, y, z));
    };
    P_2_ = [=](double t, double x, double y, double z) {
        return (alpha_2_ * b(t, x, y, z) + gamma_2_ * h(t, x, y, z));
    };
    P_3_ = [=](double t, double x, double y, double z) {
        return (alpha_3_ * c(t, x, y, z) + gamma_3_ * i(t, x, y, z));
    };
    S_1_ = [=](double t, double x, double y, double z) {
        return (alpha_1_ * a(t, x, y, z) - sixth * rho_ * j(t, x, y, z));
    };
    S_2_ = [=](double t, double x, double y, double z) {
        return (alpha_2_ * b(t, x, y, z) - sixth * rho_ * j(t, x, y, z));
    };
    S_3_ = [=](double t, double x, double y, double z) {
        return (alpha_3_ * c(t, x, y, z) - sixth * rho_ * j(t, x, y, z));
    };
    C_ = [=](double t, double x, double y, double z) { return c(t, x, y, z); };
    D_ = [=](double t, double x, double y, double z) { return d(t, x, y, z); };
    E_ = [=](double t, double x, double y, double z) { return e(t, x, y, z); };
    F_ = [=](double t, double x, double y, double z) { return f(t, x, y, z); };
    I_ = [=](double t, double x, double y, double z) { return i(t, x, y, z); };
    J_ = [=](double t, double x, double y, double z) { return j(t, x, y, z); };
    H_ = [=](double t, double x, double y, double z) { return h(t, x, y, z); };
    G_ = [=](double t, double x, double y, double z) { return g(t, x, y, z); };
}

hhw_implicit_coefficients::hhw_implicit_coefficients(heat_data_transform_3d_ptr const &heat_data_config,
                                                     pde_discretization_config_3d_ptr const &discretization_config,
                                                     splitting_method_config_ptr const splitting_config,
                                                     double const &theta)
    : theta_{theta}
{
    initialize(discretization_config, splitting_config);
    initialize_coefficients(heat_data_config);
}

} // namespace three_dimensional

} // namespace lss_pde_solvers
