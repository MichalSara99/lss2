#include "heat_coefficient_data_config_1d.hpp"

namespace lss
{

heat_coefficient_data_config_1d_builder::heat_coefficient_data_config_1d_builder()
{
}

heat_coefficient_data_config_1d_builder &heat_coefficient_data_config_1d_builder::a_coefficient(
    std::function<double(double, double)> const &a_coefficient)
{
    a_coefficient_ = a_coefficient;
    return *this;
}

heat_coefficient_data_config_1d_builder &heat_coefficient_data_config_1d_builder::b_coefficient(
    std::function<double(double, double)> const &b_coefficient)
{
    b_coefficient_ = b_coefficient;
    return *this;
}

heat_coefficient_data_config_1d_builder &heat_coefficient_data_config_1d_builder::c_coefficient(
    std::function<double(double, double)> const &c_coefficient)
{
    c_coefficient_ = c_coefficient;
    return *this;
}

heat_coefficient_data_config_1d_ptr heat_coefficient_data_config_1d_builder::build()
{
    return std::make_shared<heat_coefficient_data_config_1d>(a_coefficient_, b_coefficient_, c_coefficient_);
}

} // namespace lss
