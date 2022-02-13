#include "wave_data_config_1d.hpp"

namespace lss
{

wave_data_config_1d_builder::wave_data_config_1d_builder()
{
}

wave_data_config_1d_builder &wave_data_config_1d_builder::coefficient_data_config(
    wave_coefficient_data_config_1d_ptr const &coefficient_data_config)
{
    coefficient_data_config_ = coefficient_data_config;
    return *this;
}

wave_data_config_1d_builder &wave_data_config_1d_builder::initial_data_config(
    wave_initial_data_config_1d_ptr const &initial_data_config)
{
    initial_data_config_ = initial_data_config;
    return *this;
}

wave_data_config_1d_builder &wave_data_config_1d_builder::source_data_config(
    wave_source_data_config_1d_ptr const &source_data_config)
{
    source_data_config_ = source_data_config;
    return *this;
}

wave_data_config_1d_ptr wave_data_config_1d_builder::build()
{
    return std::make_shared<wave_data_config_1d>(coefficient_data_config_, initial_data_config_, source_data_config_);
}

} // namespace lss
