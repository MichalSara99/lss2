#include "ode_coefficient_data_config.hpp"

namespace lss
{

ode_coefficient_data_config_builder::ode_coefficient_data_config_builder()
{
}

ode_coefficient_data_config_builder &ode_coefficient_data_config_builder::a_coefficient(
    std::function<double(double)> a_coefficient)
{
    a_coeff_ = a_coefficient;
    return *this;
}

ode_coefficient_data_config_builder &ode_coefficient_data_config_builder::b_coefficient(
    std::function<double(double)> b_coefficient)
{
    b_coeff_ = b_coefficient;
    return *this;
}

ode_coefficient_data_config_ptr ode_coefficient_data_config_builder::build()
{
    return std::make_shared<ode_coefficient_data_config>(a_coeff_, b_coeff_);
}

} // namespace lss
