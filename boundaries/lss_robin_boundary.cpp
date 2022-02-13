#include "lss_robin_boundary.hpp"

namespace lss_boundary
{

robin_boundary_1d::robin_boundary_1d(double linear_value, double value) : boundary_1d(linear_value, value)
{
}

robin_boundary_1d::robin_boundary_1d(const std::function<double(double)> &linear_value,
                                     const std::function<double(double)> &value)
    : boundary_1d(linear_value, value)
{
}

double robin_boundary_1d::linear_value() const
{
    return this->linear_val_;
}
double robin_boundary_1d::value() const
{
    return this->const_val_;
}

double robin_boundary_1d::linear_value(double time) const
{
    LSS_ASSERT(is_time_dependent_ == true, "robin_boundary_1d: Boundary must not be time independent.");
    return this->linear_fun_(time);
}
double robin_boundary_1d::value(double time) const
{
    LSS_ASSERT(is_time_dependent_ == true, "robin_boundary_1d: Boundary must not be time independent.");
    return this->const_fun_(time);
}

robin_boundary_2d::robin_boundary_2d(const std::function<double(double, double)> &linear_value,
                                     const std::function<double(double, double)> &value)
    : boundary_2d(linear_value, value)
{
}

double robin_boundary_2d::linear_value(double time, double space_arg) const
{
    return this->linear_(time, space_arg);
}
double robin_boundary_2d::value(double time, double space_arg) const
{
    return this->const_(time, space_arg);
}

robin_boundary_3d::robin_boundary_3d(const std::function<double(double, double, double)> &linear_value,
                                     const std::function<double(double, double, double)> &value)
    : boundary_3d(linear_value, value)
{
}

double robin_boundary_3d::linear_value(double time, double space_1_arg, double space_2_arg) const
{
    return this->linear_(time, space_1_arg, space_2_arg);
}
double robin_boundary_3d::value(double time, double space_1_arg, double space_2_arg) const
{
    return this->const_(time, space_1_arg, space_2_arg);
}
} // namespace lss_boundary
