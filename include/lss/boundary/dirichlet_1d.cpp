#include "dirichlet_1d.hpp"

namespace lss
{

dirichlet_1d_builder::dirichlet_1d_builder()
{
}

dirichlet_1d_builder &dirichlet_1d_builder::value(const std::function<double(double)> &value)
{
    fun_value_ = value;
    return *this;
}

dirichlet_1d_builder &dirichlet_1d_builder::value(double value)
{
    const_value_ = value;
    return *this;
}

dirichlet_1d_ptr dirichlet_1d_builder::build()
{
    if (fun_value_ != nullptr)
        return std::make_shared<dirichlet_1d>(fun_value_);
    else
        return std::make_shared<dirichlet_1d>(const_value_);
}

} // namespace lss
