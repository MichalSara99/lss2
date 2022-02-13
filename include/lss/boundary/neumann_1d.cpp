#include "neumann_1d.hpp"

namespace lss
{

neumann_1d_builder::neumann_1d_builder()
{
}

neumann_1d_builder &neumann_1d_builder::value(const std::function<double(double)> &value)
{
    fun_value_ = value;
    return *this;
}

neumann_1d_builder &neumann_1d_builder::value(double value)
{
    const_value_ = value;
    return *this;
}

neumann_1d_ptr neumann_1d_builder::build()
{
    if (fun_value_ != nullptr)
        return std::make_shared<neumann_1d>(fun_value_);
    else
        return std::make_shared<neumann_1d>(const_value_);
}

} // namespace lss
