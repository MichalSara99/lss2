#if !defined(_IMPLICIT_PDE_SCHEME_HPP_)
#define _IMPLICIT_PDE_SCHEME_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_implicit_pde_scheme.hpp"

namespace lss
{

using implicit_pde_schemes = lss_enumerations::implicit_pde_schemes_enum;
using implicit_pde_scheme_ptr = lss_pde_solvers::implicit_pde_scheme_ptr;
using implicit_pde_scheme = lss_pde_solvers::implicit_pde_scheme;

struct implicit_pde_scheme_builder
{
  private:
    double value_ = -1.0;
    implicit_pde_schemes scheme_;

  public:
    LSS_API explicit implicit_pde_scheme_builder();

    LSS_API implicit_pde_scheme_builder &value(double value);

    LSS_API implicit_pde_scheme_builder &value(implicit_pde_schemes schemes_enum);

    LSS_API implicit_pde_scheme_ptr build();
};

} // namespace lss

#endif ///_IMPLICIT_PDE_SCHEME_HPP_
