#include "lss_implicit_pde_scheme.hpp"

namespace lss_pde_solvers
{

implicit_pde_scheme::implicit_pde_scheme()
{
}

implicit_pde_scheme::implicit_pde_scheme(double value) : value_{value}
{
    LSS_ASSERT(((value_ >= 0.0) && (value_ <= 1.0)), "implicit_pde_scheme: value must be in <0.0,1.0>.");
}

implicit_pde_scheme::implicit_pde_scheme(implicit_pde_schemes_enum implicit_pde_enum)
    : value_{(implicit_pde_enum == implicit_pde_schemes_enum::CrankNicolson) ? 0.5 : 1.0}
{
}

double implicit_pde_scheme::value() const
{
    return value_;
}

implicit_pde_scheme::~implicit_pde_scheme()
{
}

} // namespace lss_pde_solvers
