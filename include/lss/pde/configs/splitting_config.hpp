#if !defined(_SPLITTING_CONFIG_HPP_)
#define _SPLITTING_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_splitting_method_config.hpp"

namespace lss
{
using splitting_config_ptr = lss_pde_solvers::splitting_method_config_ptr;
using splitting_config = lss_pde_solvers::splitting_method_config;
using splitting_method = lss_enumerations::splitting_method_enum;

struct splitting_config_builder
{
  private:
    splitting_method splitting_method_;
    double weighting_value_;

  public:
    LSS_API explicit splitting_config_builder();

    LSS_API splitting_config_builder &splitting(splitting_method splitting_method);

    LSS_API splitting_config_builder &weighting_value(double value);

    LSS_API splitting_config_ptr build();
};

} // namespace lss

#endif ///_SPLITTING_CONFIG_HPP_
