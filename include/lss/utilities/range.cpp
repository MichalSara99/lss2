#include "range.hpp"

namespace lss
{

range_builder::range_builder()
{
}

range_builder &range_builder::lower(double value)
{
    l_ = value;
    return *this;
}

range_builder &range_builder::upper(double value)
{
    u_ = value;
    return *this;
}

lss_utility::range_ptr range_builder::build()
{
    return std::make_shared<lss_utility::range>(l_, u_);
}

} // namespace lss
