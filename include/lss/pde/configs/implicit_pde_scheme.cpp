#include "implicit_pde_scheme.hpp"

namespace lss
{

implicit_pde_scheme_builder::implicit_pde_scheme_builder()
{
}

implicit_pde_scheme_builder &implicit_pde_scheme_builder::value(double value)
{
    value_ = value;
    return *this;
}

implicit_pde_scheme_builder &implicit_pde_scheme_builder::value(implicit_pde_schemes schemes_enum)
{
    scheme_ = schemes_enum;
    return *this;
}

implicit_pde_scheme_ptr implicit_pde_scheme_builder::build()
{
    // prioritise the explicit value:
    if (value_ != -1.0)
        return std::make_shared<implicit_pde_scheme>(value_);
    else
        return std::make_shared<implicit_pde_scheme>(scheme_);
}
} // namespace lss
