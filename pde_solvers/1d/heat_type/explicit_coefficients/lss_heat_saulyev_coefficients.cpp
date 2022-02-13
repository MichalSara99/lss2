#include "lss_heat_saulyev_coefficients.hpp"

namespace lss_pde_solvers
{

namespace one_dimensional
{

void heat_saulyev_coefficients::initialize_coefficients(heat_coefficients_ptr const &coefficients)
{
    const double one = 1.0;
    // save coefficients locally:
    auto const a = coefficients->A_;
    auto const b = coefficients->B_;
    auto const d = coefficients->D_;
    k_ = coefficients->k_;

    A_ = [=](double t, double x) { return (a(t, x) / (one + b(t, x))); };
    B_ = [=](double t, double x) { return ((one - b(t, x)) / (one + b(t, x))); };
    D_ = [=](double t, double x) { return (d(t, x) / (one + b(t, x))); };
    K_ = [=](double t, double x) { return (k_ / (one + b(t, x))); };

    space_size_ = coefficients->space_size_;
}

heat_saulyev_coefficients::heat_saulyev_coefficients(heat_coefficients_ptr const &coefficients)
{
    initialize_coefficients(coefficients);
}

} // namespace one_dimensional

} // namespace lss_pde_solvers
