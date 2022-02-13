#include "heat_initial_data_config_2d.hpp"

namespace lss
{

heat_initial_data_config_2d_builder::heat_initial_data_config_2d_builder()
{
}

heat_initial_data_config_2d_builder &heat_initial_data_config_2d_builder::condition(
    std::function<double(double, double)> const &condition)
{
    initial_condition_ = condition;
    return *this;
}

heat_initial_data_config_2d_ptr heat_initial_data_config_2d_builder::build()
{
    return std::make_shared<heat_initial_data_config_2d>(initial_condition_);
}

} // namespace lss
