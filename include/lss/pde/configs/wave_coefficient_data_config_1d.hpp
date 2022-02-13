#if !defined(_WAVE_COEFFICIENT_DATA_CONFIG_1D_HPP_)
#define _WAVE_COEFFICIENT_DATA_CONFIG_1D_HPP_

#include <functional>

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_wave_data_config.hpp"

namespace lss
{

using wave_coefficient_data_config_1d = lss_pde_solvers::wave_coefficient_data_config_1d;

using wave_coefficient_data_config_1d_ptr = lss_pde_solvers::sptr_t<lss_pde_solvers::wave_coefficient_data_config_1d>;

struct wave_coefficient_data_config_1d_builder
{
  private:
    std::function<double(double, double)> a_coefficient_;
    std::function<double(double, double)> b_coefficient_;
    std::function<double(double, double)> c_coefficient_;
    std::function<double(double, double)> d_coefficient_;

  public:
    LSS_API explicit wave_coefficient_data_config_1d_builder();

    LSS_API wave_coefficient_data_config_1d_builder &a_coefficient(
        std::function<double(double, double)> const &a_coefficient);

    LSS_API wave_coefficient_data_config_1d_builder &b_coefficient(
        std::function<double(double, double)> const &b_coefficient);

    LSS_API wave_coefficient_data_config_1d_builder &c_coefficient(
        std::function<double(double, double)> const &c_coefficient);

    LSS_API wave_coefficient_data_config_1d_builder &d_coefficient(
        std::function<double(double, double)> const &d_coefficient);

    LSS_API wave_coefficient_data_config_1d_ptr build();
};

} // namespace lss

#endif ///_WAVE_COEFFICIENT_DATA_CONFIG_1D_HPP_
