#include "heat_initial_data_config_1d.hpp"

namespace lss
{

heat_initial_data_config_1d_builder::heat_initial_data_config_1d_builder()
{
}

heat_initial_data_config_1d_builder &heat_initial_data_config_1d_builder::condition(
    std::function<double(double)> const &condition)
{
    initial_condition_ = condition;
    return *this;
}

heat_initial_data_config_1d_ptr heat_initial_data_config_1d_builder::build()
{
    return std::make_shared<heat_initial_data_config_1d>(initial_condition_);
}

} // namespace lss
