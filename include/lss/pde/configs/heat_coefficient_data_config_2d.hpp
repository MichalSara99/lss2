#if !defined(_HEAT_COEFFICIENT_DATA_CONFIG_2D_HPP_)
#define _HEAT_COEFFICIENT_DATA_CONFIG_2D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_data_config.hpp"

namespace lss
{

using heat_coefficient_data_config_2d = lss_pde_solvers::heat_coefficient_data_config_2d;
using heat_coefficient_data_config_2d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::heat_coefficient_data_config_2d>;

struct heat_coefficient_data_config_2d_builder
{
  private:
    std::function<double(double, double, double)> a_coefficient_;
    std::function<double(double, double, double)> b_coefficient_;
    std::function<double(double, double, double)> c_coefficient_;
    std::function<double(double, double, double)> d_coefficient_;
    std::function<double(double, double, double)> e_coefficient_;
    std::function<double(double, double, double)> f_coefficient_;

  public:
    LSS_API explicit heat_coefficient_data_config_2d_builder();

    LSS_API heat_coefficient_data_config_2d_builder &a_coefficient(
        std::function<double(double, double, double)> const &a_coefficient);

    LSS_API heat_coefficient_data_config_2d_builder &b_coefficient(
        std::function<double(double, double, double)> const &b_coefficient);

    LSS_API heat_coefficient_data_config_2d_builder &c_coefficient(
        std::function<double(double, double, double)> const &c_coefficient);

    LSS_API heat_coefficient_data_config_2d_builder &d_coefficient(
        std::function<double(double, double, double)> const &d_coefficient);

    LSS_API heat_coefficient_data_config_2d_builder &e_coefficient(
        std::function<double(double, double, double)> const &e_coefficient);

    LSS_API heat_coefficient_data_config_2d_builder &f_coefficient(
        std::function<double(double, double, double)> const &f_coefficient);

    LSS_API heat_coefficient_data_config_2d_ptr build();
};

} // namespace lss

#endif ///_HEAT_COEFFICIENT_DATA_CONFIG_2D_HPP_
