#include "lss_heat_euler_coefficients_2d.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

void heat_euler_coefficients_2d::initialize_coefficients(heat_coefficients_2d_ptr const &coefficients)
{
    // time step:
    k_ = coefficients->k_;
    // size of spaces discretization:
    space_size_x_ = coefficients->space_size_x_;
    space_size_y_ = coefficients->space_size_y_;
    // calculate scheme coefficients:
    rho_ = coefficients->rho_;
    // save coefficients locally:
    M_ = coefficients->M_;
    M_tilde_ = coefficients->M_tilde_;
    P_ = coefficients->P_;
    P_tilde_ = coefficients->P_tilde_;
    Z_ = coefficients->Z_;
    W_ = coefficients->W_;
    auto const &c = coefficients->C_;
    auto const gamma = coefficients->gamma_;
    C_ = [=](double t, double x, double y) { return (gamma * c(t, x, y)); };
}

heat_euler_coefficients_2d::heat_euler_coefficients_2d(heat_coefficients_2d_ptr const &coefficients)
{
    initialize_coefficients(coefficients);
}

} // namespace two_dimensional

} // namespace lss_pde_solvers
