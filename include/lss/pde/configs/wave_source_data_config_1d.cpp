#include "wave_source_data_config_1d.hpp"

namespace lss
{

wave_source_data_config_1d_builder::wave_source_data_config_1d_builder()
{
}

wave_source_data_config_1d_builder &wave_source_data_config_1d_builder::wave_source(
    std::function<double(double, double)> const &wave_source)
{
    wave_source_ = wave_source;
    return *this;
}

wave_source_data_config_1d_ptr wave_source_data_config_1d_builder::build()
{
    return std::make_shared<wave_source_data_config_1d>(wave_source_);
}

} // namespace lss
