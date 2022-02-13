#include "lss_boundary.hpp"

namespace lss_boundary
{
boundary_1d::boundary_1d(double linear, double constant)
    : linear_val_{linear}, const_val_{constant}, is_time_dependent_{false}
{
}

boundary_1d::boundary_1d(const std::function<double(double)> &linear, const std::function<double(double)> &constant)
    : linear_fun_{linear}, const_fun_{constant}, is_time_dependent_{true}
{
}

boundary_1d::~boundary_1d()
{
}

bool const &boundary_1d::is_time_dependent() const
{
    return is_time_dependent_;
}

boundary_2d::boundary_2d(const std::function<double(double, double)> &linear,
                         const std::function<double(double, double)> &constant)
    : linear_{linear}, const_{constant}
{
}

boundary_2d ::~boundary_2d()
{
}

boundary_3d::boundary_3d(const std::function<double(double, double, double)> &linear,
                         const std::function<double(double, double, double)> &constant)
    : linear_{linear}, const_{constant}
{
}

boundary_3d ::~boundary_3d()
{
}

} // namespace lss_boundary
