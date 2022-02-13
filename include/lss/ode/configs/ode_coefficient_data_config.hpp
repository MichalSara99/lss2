#if !defined(_ODE_COEFFICIENT_DATA_CONFIG_HPP_)
#define _ODE_COEFFICIENT_DATA_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../ode_solvers/lss_ode_data_config.hpp"

namespace lss
{

using ode_coefficient_data_config_ptr = lss_ode_solvers::ode_coefficient_data_config_ptr;
using ode_coefficient_data_config = lss_ode_solvers::ode_coefficient_data_config;

struct ode_coefficient_data_config_builder
{
  private:
    std::function<double(double)> a_coeff_;
    std::function<double(double)> b_coeff_;

  public:
    LSS_API explicit ode_coefficient_data_config_builder();

    LSS_API ode_coefficient_data_config_builder &a_coefficient(std::function<double(double)> a_coefficient);

    LSS_API ode_coefficient_data_config_builder &b_coefficient(std::function<double(double)> b_coefficient);

    LSS_API ode_coefficient_data_config_ptr build();
};

} // namespace lss

#endif ///_ODE_COEFFICIENT_DATA_CONFIG_HPP_
