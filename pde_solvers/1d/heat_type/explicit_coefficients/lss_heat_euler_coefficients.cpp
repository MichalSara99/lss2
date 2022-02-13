#include "lss_heat_euler_coefficients.hpp"
namespace lss_pde_solvers
{

namespace one_dimensional
{

void heat_euler_coefficients::initialize_coefficients(heat_coefficients_ptr const &coefficients)
{
    // save coefficients locally:
    auto const a = coefficients->A_;
    auto const b = coefficients->B_;
    auto const d = coefficients->D_;

    const double one = 1.0;
    const double two = 2.0;

    A_ = [=](double t, double x) { return a(t, x); };
    B_ = [=](double t, double x) { return (one - two * b(t, x)); };
    D_ = [=](double t, double x) { return d(t, x); };
    k_ = coefficients->k_;
    space_size_ = coefficients->space_size_;
}

heat_euler_coefficients::heat_euler_coefficients(heat_coefficients_ptr const &coefficients)
{
    initialize_coefficients(coefficients);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
