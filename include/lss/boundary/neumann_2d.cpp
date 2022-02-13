#include "neumann_2d.hpp"

namespace lss
{

neumann_2d_builder::neumann_2d_builder()
{
}

neumann_2d_builder &neumann_2d_builder::value(const std::function<double(double, double)> &value)
{
    value_ = value;
    return *this;
}

neumann_2d_ptr neumann_2d_builder::build()
{
    return std::make_shared<neumann_2d>(value_);
}

} // namespace lss
