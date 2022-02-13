#include "heat_coefficient_data_config_2d.hpp"

namespace lss
{

heat_coefficient_data_config_2d_builder::heat_coefficient_data_config_2d_builder()
{
}

heat_coefficient_data_config_2d_builder &heat_coefficient_data_config_2d_builder::a_coefficient(
    std::function<double(double, double, double)> const &a_coefficient)
{
    a_coefficient_ = a_coefficient;
    return *this;
}

heat_coefficient_data_config_2d_builder &heat_coefficient_data_config_2d_builder::b_coefficient(
    std::function<double(double, double, double)> const &b_coefficient)
{
    b_coefficient_ = b_coefficient;
    return *this;
}

heat_coefficient_data_config_2d_builder &heat_coefficient_data_config_2d_builder::c_coefficient(
    std::function<double(double, double, double)> const &c_coefficient)
{
    c_coefficient_ = c_coefficient;
    return *this;
}

heat_coefficient_data_config_2d_builder &heat_coefficient_data_config_2d_builder::d_coefficient(
    std::function<double(double, double, double)> const &d_coefficient)
{
    d_coefficient_ = d_coefficient;
    return *this;
}

heat_coefficient_data_config_2d_builder &heat_coefficient_data_config_2d_builder::e_coefficient(
    std::function<double(double, double, double)> const &e_coefficient)
{
    e_coefficient_ = e_coefficient;
    return *this;
}

heat_coefficient_data_config_2d_builder &heat_coefficient_data_config_2d_builder::f_coefficient(
    std::function<double(double, double, double)> const &f_coefficient)
{
    f_coefficient_ = f_coefficient;
    return *this;
}

heat_coefficient_data_config_2d_ptr heat_coefficient_data_config_2d_builder::build()
{
    return std::make_shared<heat_coefficient_data_config_2d>(a_coefficient_, b_coefficient_, c_coefficient_,
                                                             d_coefficient_, e_coefficient_, f_coefficient_);
}

} // namespace lss
