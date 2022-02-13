#include "ode_data_config.hpp"

namespace lss
{

ode_data_config_builder::ode_data_config_builder()
{
}

ode_data_config_builder &ode_data_config_builder::coefficient_data_config(
    ode_coefficient_data_config_ptr const &coefficient_data_config)
{
    coefficient_data_cfg_ = coefficient_data_config;
    return *this;
}

ode_data_config_builder &ode_data_config_builder::nonhom_data_config(
    ode_nonhom_data_config_ptr const &nonhom_data_config)
{
    nonhom_data_cfg_ = nonhom_data_config;
    return *this;
}

ode_data_config_ptr ode_data_config_builder::build()
{
    return std::make_shared<ode_data_config>(coefficient_data_cfg_, nonhom_data_cfg_);
}

} // namespace lss
