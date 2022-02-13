#include "heat_data_config_1d.hpp"

namespace lss
{

heat_data_config_1d_builder::heat_data_config_1d_builder()
{
}

heat_data_config_1d_builder &heat_data_config_1d_builder::coefficient_data_config(
    heat_coefficient_data_config_1d_ptr const &coefficient_data_config)
{
    coefficient_data_config_ = coefficient_data_config;
    return *this;
}

heat_data_config_1d_builder &heat_data_config_1d_builder::initial_data_config(
    heat_initial_data_config_1d_ptr const &initial_data_config)
{
    initial_data_config_ = initial_data_config;
    return *this;
}

heat_data_config_1d_builder &heat_data_config_1d_builder::source_data_config(
    heat_source_data_config_1d_ptr const &source_data_config)
{
    source_data_config_ = source_data_config;
    return *this;
}

heat_data_config_1d_ptr heat_data_config_1d_builder::build()
{
    return std::make_shared<heat_data_config_1d>(coefficient_data_config_, initial_data_config_, source_data_config_);
}

} // namespace lss
