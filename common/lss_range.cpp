#include "lss_range.hpp"

namespace lss_utility
{

range::range(double lower, double upper) : l_{lower}, u_{upper}
{
}

range::range() : range(double{}, double{})
{
}

range::~range()
{
}

range::range(range const &copy) : l_{copy.l_}, u_{copy.u_}
{
}
range::range(range &&other) noexcept : l_{std::move(other.l_)}, u_{std::move(other.u_)}
{
}

range &range::operator=(range const &copy)
{
    if (this != &copy)
    {
        l_ = copy.l_;
        u_ = copy.u_;
    }
    return *this;
}

range &range::operator=(range &&other) noexcept
{
    if (this != &other)
    {
        l_ = std::move(other.l_);
        u_ = std::move(other.u_);
    }
    return *this;
}

double range::lower() const
{
    return l_;
}
double range::upper() const
{
    return u_;
}

double range::spread() const
{
    return (u_ - l_);
}

double range::mid_point() const
{
    return 0.5 * (l_ + u_);
}

} // namespace lss_utility
