#include "wave_coefficient_data_config_1d.hpp"

namespace lss
{

wave_coefficient_data_config_1d_builder::wave_coefficient_data_config_1d_builder()
{
}

wave_coefficient_data_config_1d_builder &wave_coefficient_data_config_1d_builder::a_coefficient(
    std::function<double(double, double)> const &a_coefficient)
{
    a_coefficient_ = a_coefficient;
    return *this;
}

wave_coefficient_data_config_1d_builder &wave_coefficient_data_config_1d_builder::b_coefficient(
    std::function<double(double, double)> const &b_coefficient)
{
    b_coefficient_ = b_coefficient;
    return *this;
}

wave_coefficient_data_config_1d_builder &wave_coefficient_data_config_1d_builder::c_coefficient(
    std::function<double(double, double)> const &c_coefficient)
{
    c_coefficient_ = c_coefficient;
    return *this;
}

wave_coefficient_data_config_1d_builder &wave_coefficient_data_config_1d_builder::d_coefficient(
    std::function<double(double, double)> const &d_coefficient)
{
    d_coefficient_ = d_coefficient;
    return *this;
}

wave_coefficient_data_config_1d_ptr wave_coefficient_data_config_1d_builder::build()
{
    return std::make_shared<wave_coefficient_data_config_1d>(a_coefficient_, b_coefficient_, c_coefficient_,
                                                             d_coefficient_);
}

} // namespace lss
