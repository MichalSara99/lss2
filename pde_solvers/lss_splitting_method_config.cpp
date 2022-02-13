#include "lss_splitting_method_config.hpp"

namespace lss_pde_solvers
{

splitting_method_config::splitting_method_config(splitting_method_enum splittitng_method, double weighting_value)
    : splitting_method_{splittitng_method}, weighting_value_{weighting_value}
{
}

splitting_method_config::~splitting_method_config()
{
}

splitting_method_enum splitting_method_config::splitting_method() const
{
    return splitting_method_;
}

double splitting_method_config::weighting_value() const
{
    return weighting_value_;
}

} // namespace lss_pde_solvers
