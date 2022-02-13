#if !defined(_ODE_DATA_CONFIG_HPP_)
#define _ODE_DATA_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../ode_solvers/lss_ode_data_config.hpp"
#include "../../../lss/ode/configs/ode_coefficient_data_config.hpp"
#include "../../../lss/ode/configs/ode_nonhom_data_config.hpp"

namespace lss
{

using ode_data_config_ptr = lss_ode_solvers::ode_data_config_ptr;
using ode_data_config = lss_ode_solvers::ode_data_config;

struct ode_data_config_builder
{
  private:
    ode_coefficient_data_config_ptr coefficient_data_cfg_;
    ode_nonhom_data_config_ptr nonhom_data_cfg_;

  public:
    LSS_API explicit ode_data_config_builder();

    LSS_API ode_data_config_builder &coefficient_data_config(
        ode_coefficient_data_config_ptr const &coefficient_data_config);

    LSS_API ode_data_config_builder &nonhom_data_config(ode_nonhom_data_config_ptr const &nonhom_data_config);

    LSS_API ode_data_config_ptr build();
};

} // namespace lss

#endif ///_ODE_DATA_CONFIG_HPP_
