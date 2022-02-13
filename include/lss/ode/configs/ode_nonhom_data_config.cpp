#include "ode_nonhom_data_config.hpp"

namespace lss
{

ode_nonhom_data_config_builder::ode_nonhom_data_config_builder()
{
}

ode_nonhom_data_config_builder &ode_nonhom_data_config_builder::nonhom_function(
    std::function<double(double)> const &function)
{
    nonhom_fun_ = function;
    return *this;
}

ode_nonhom_data_config_ptr ode_nonhom_data_config_builder::build()
{
    return std::make_shared<ode_nonhom_data_config>(nonhom_fun_);
}

} // namespace lss
