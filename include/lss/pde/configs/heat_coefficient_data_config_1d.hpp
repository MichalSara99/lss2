#if !defined(_HEAT_COEFFICIENT_DATA_CONFIG_1D_HPP_)
#define _HEAT_COEFFICIENT_DATA_CONFIG_1D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_data_config.hpp"

namespace lss
{

using heat_coefficient_data_config_1d = lss_pde_solvers::heat_coefficient_data_config_1d;

using heat_coefficient_data_config_1d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::heat_coefficient_data_config_1d>;

struct heat_coefficient_data_config_1d_builder
{
  private:
    std::function<double(double, double)> a_coefficient_;
    std::function<double(double, double)> b_coefficient_;
    std::function<double(double, double)> c_coefficient_;

  public:
    LSS_API explicit heat_coefficient_data_config_1d_builder();

    LSS_API heat_coefficient_data_config_1d_builder &a_coefficient(
        std::function<double(double, double)> const &a_coefficient);

    LSS_API heat_coefficient_data_config_1d_builder &b_coefficient(
        std::function<double(double, double)> const &b_coefficient);

    LSS_API heat_coefficient_data_config_1d_builder &c_coefficient(
        std::function<double(double, double)> const &c_coefficient);

    LSS_API heat_coefficient_data_config_1d_ptr build();
};

} // namespace lss

#endif ///_HEAT_COEFFICIENT_DATA_CONFIG_1D_HPP_
