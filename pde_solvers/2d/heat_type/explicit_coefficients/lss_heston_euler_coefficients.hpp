#if !defined(_LSS_HESTON_EULER_COEFFICIENTS_HPP_)
#define _LSS_HESTON_EULER_COEFFICIENTS_HPP_

#include "../../../../common/lss_utility.hpp"
#include "../implicit_coefficients/lss_heston_implicit_coefficients.hpp"

namespace lss_pde_solvers
{
namespace two_dimensional
{

using lss_utility::sptr_t;

/**
    heston_euler_coefficients object
 */
struct heston_euler_coefficients
{
  public:
    // scheme constant coefficients:
    double rho_, k_;
    std::size_t space_size_x_, space_size_y_;
    // functional coefficients:
    std::function<double(double, double, double)> M_;
    std::function<double(double, double, double)> M_tilde_;
    std::function<double(double, double, double)> P_;
    std::function<double(double, double, double)> P_tilde_;
    std::function<double(double, double, double)> Z_;
    std::function<double(double, double, double)> W_;
    std::function<double(double, double, double)> C_;

  private:
    void initialize_coefficients(heston_implicit_coefficients_ptr const &coefficients);

  public:
    heston_euler_coefficients() = delete;

    heston_euler_coefficients(heston_implicit_coefficients_ptr const &coefficients);
};

using heston_euler_coefficients_ptr = sptr_t<heston_euler_coefficients>;

} // namespace two_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HESTON_EULER_COEFFICIENTS_HPP_
