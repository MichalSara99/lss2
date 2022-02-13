#include "robin_1d.hpp"

namespace lss
{

robin_1d_builder::robin_1d_builder()
{
}

robin_1d_builder &robin_1d_builder::values(const std::function<double(double)> &linear_value,
                                           const std::function<double(double)> &value)
{
    fun_linear_value_ = linear_value;
    fun_value_ = value;
    return *this;
}

robin_1d_builder &robin_1d_builder::values(double linear_value, double value)
{
    const_linear_value_ = linear_value;
    const_value_ = value;
    return *this;
}

robin_1d_ptr robin_1d_builder::build()
{
    if ((fun_value_ != nullptr) && (fun_linear_value_ != nullptr))
    {
        return std::make_shared<robin_1d>(fun_linear_value_, fun_value_);
    }
    else
    {
        return std::make_shared<robin_1d>(const_linear_value_, const_value_);
    }
}

} // namespace lss
