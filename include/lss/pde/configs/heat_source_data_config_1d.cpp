#include "heat_source_data_config_1d.hpp"

namespace lss
{

heat_source_data_config_1d_builder::heat_source_data_config_1d_builder()
{
}

heat_source_data_config_1d_builder &heat_source_data_config_1d_builder::heat_source(
    std::function<double(double, double)> const &heat_source)
{
    heat_source_ = heat_source;
    return *this;
}

heat_source_data_config_1d_ptr heat_source_data_config_1d_builder::build()
{
    return std::make_shared<heat_source_data_config_1d>(heat_source_);
}

} // namespace lss
