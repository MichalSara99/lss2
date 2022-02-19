#include "dirichlet_2d.hpp"

namespace lss
{

dirichlet_2d_builder::dirichlet_2d_builder()
{
}

dirichlet_2d_builder &dirichlet_2d_builder::value(const std::function<double(double, double)> &value)
{
    value_ = value;
    return *this;
}

dirichlet_2d_ptr dirichlet_2d_builder::build()
{
    return std::make_shared<dirichlet_2d>(value_);
}

} // namespace lss
