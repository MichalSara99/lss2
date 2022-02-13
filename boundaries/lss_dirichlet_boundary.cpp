#include "lss_dirichlet_boundary.hpp"

namespace lss_boundary
{
dirichlet_boundary_1d::dirichlet_boundary_1d(const std::function<double(double)> &value) : boundary_1d(nullptr, value)
{
}

dirichlet_boundary_1d::dirichlet_boundary_1d(double value) : boundary_1d(double{}, value)
{
}

double dirichlet_boundary_1d::value(double time) const
{
    LSS_ASSERT(is_time_dependent_ == true, "dirichlet_boundary_1d: Boundary must not be time independent.");
    return this->const_fun_(time);
}

double dirichlet_boundary_1d::value() const
{
    return this->const_val_;
}

dirichlet_boundary_2d::dirichlet_boundary_2d(const std::function<double(double, double)> &value)
    : boundary_2d(nullptr, value)
{
}

double dirichlet_boundary_2d::value(double time, double space_arg) const
{
    return this->const_(time, space_arg);
}

dirichlet_boundary_3d::dirichlet_boundary_3d(const std::function<double(double, double, double)> &value)
    : boundary_3d(nullptr, value)
{
}

double dirichlet_boundary_3d::value(double time, double space_1_arg, double space_2_arg) const
{
    return this->const_(time, space_1_arg, space_2_arg);
}

} // namespace lss_boundary
