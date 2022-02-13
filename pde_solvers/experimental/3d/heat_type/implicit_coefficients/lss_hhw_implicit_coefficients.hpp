#if !defined(_LSS_HHW_IMPLICIT_COEFFICIENTS_HPP_)
#define _LSS_HHW_IMPLICIT_COEFFICIENTS_HPP_

#include "../../../../../common/lss_utility.hpp"
#include "../../../../../discretization/lss_discretization.hpp"
#include "../../../../lss_pde_discretization_config.hpp"
#include "../../../../lss_splitting_method_config.hpp"
#include "../../../../transformation/lss_heat_data_transform.hpp"

namespace lss_pde_solvers
{
namespace three_dimensional
{

using lss_utility::range_ptr;
using lss_utility::sptr_t;

/**
    hhw_implicit_coefficients object
 */
struct hhw_implicit_coefficients
{
  public:
    // scheme constant coefficients:
    double alpha_1_, alpha_2_, alpha_3_;
    double beta_1_, beta_2_, beta_3_;
    double gamma_1_, gamma_2_, gamma_3_;
    double k_, rho_, zeta_;
    std::size_t space_size_x_, space_size_y_, space_size_z_;
    range_ptr rangex_, rangey_, rangez_;
    // theta variable:
    double theta_;
    // functional coefficients:
    std::function<double(double, double, double, double)> M_1_;
    std::function<double(double, double, double, double)> M_2_;
    std::function<double(double, double, double, double)> M_3_;
    std::function<double(double, double, double, double)> P_1_;
    std::function<double(double, double, double, double)> P_2_;
    std::function<double(double, double, double, double)> P_3_;
    std::function<double(double, double, double, double)> S_1_;
    std::function<double(double, double, double, double)> S_2_;
    std::function<double(double, double, double, double)> S_3_;
    std::function<double(double, double, double, double)> C_;
    std::function<double(double, double, double, double)> D_;
    std::function<double(double, double, double, double)> E_;
    std::function<double(double, double, double, double)> F_;
    std::function<double(double, double, double, double)> I_;
    std::function<double(double, double, double, double)> H_;
    std::function<double(double, double, double, double)> J_;
    std::function<double(double, double, double, double)> G_;

  private:
    void initialize(pde_discretization_config_3d_ptr const &discretization_config,
                    splitting_method_config_ptr const &splitting_config);

    void initialize_coefficients(heat_data_transform_3d_ptr const &heat_data_config);

  public:
    hhw_implicit_coefficients() = delete;

    hhw_implicit_coefficients(heat_data_transform_3d_ptr const &heat_data_config,
                              pde_discretization_config_3d_ptr const &discretization_config,
                              splitting_method_config_ptr const splitting_config, double const &theta);
};

using hhw_implicit_coefficients_ptr = sptr_t<hhw_implicit_coefficients>;

} // namespace three_dimensional

} // namespace lss_pde_solvers

#endif ///_LSS_HHW_IMPLICIT_COEFFICIENTS_HPP_
