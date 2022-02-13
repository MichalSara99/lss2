#include "robin_2d.hpp"

namespace lss
{

robin_2d_builder::robin_2d_builder()
{
}

robin_2d_builder &robin_2d_builder::value(const std::function<double(double, double)> &value)
{
    value_ = value;
    return *this;
}

robin_2d_builder &robin_2d_builder::linear_value(const std::function<double(double, double)> &linear_value)
{
    linear_value_ = linear_value;
    return *this;
}

robin_2d_ptr robin_2d_builder::build()
{
    return std::make_shared<robin_2d>(linear_value_, value_);
}

} // namespace lss
