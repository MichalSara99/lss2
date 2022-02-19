#include "neumann_3d.hpp"

namespace lss
{

neumann_3d_builder::neumann_3d_builder()
{
}

neumann_3d_builder &neumann_3d_builder::value(const std::function<double(double, double, double)> &value)
{
    value_ = value;
    return *this;
}

neumann_3d_ptr neumann_3d_builder::build()
{
    return std::make_shared<neumann_3d>(value_);
}

} // namespace lss
