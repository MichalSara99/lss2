#include "robin_3d.hpp"

namespace lss
{

robin_3d_builder::robin_3d_builder()
{
}

robin_3d_builder &robin_3d_builder::value(const std::function<double(double, double, double)> &value)
{
    value_ = value;
    return *this;
}

robin_3d_builder &robin_3d_builder::linear_value(const std::function<double(double, double, double)> &linear_value)
{
    linear_value_ = linear_value;
    return *this;
}

robin_3d_ptr robin_3d_builder::build()
{
    return std::make_shared<robin_3d>(linear_value_, value_);
}

} // namespace lss
