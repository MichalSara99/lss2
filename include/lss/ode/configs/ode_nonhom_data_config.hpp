#if !defined(_ODE_NONHOM_DATA_CONFIG_HPP_)
#define _ODE_NONHOM_DATA_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../ode_solvers/lss_ode_data_config.hpp"

namespace lss
{

using ode_nonhom_data_config_ptr = lss_ode_solvers::ode_nonhom_data_config_ptr;
using ode_nonhom_data_config = lss_ode_solvers::ode_nonhom_data_config;

struct ode_nonhom_data_config_builder
{
  private:
    std::function<double(double)> nonhom_fun_;

  public:
    LSS_API explicit ode_nonhom_data_config_builder();

    LSS_API ode_nonhom_data_config_builder &nonhom_function(std::function<double(double)> const &function);

    LSS_API ode_nonhom_data_config_ptr build();
};

} // namespace lss

#endif ///_ODE_NONHOM_DATA_CONFIG_HPP_
