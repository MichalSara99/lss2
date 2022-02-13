#if !defined(_HEAT_EXPLICIT_SOLVER_CONFIG_HPP_)
#define _HEAT_EXPLICIT_SOLVER_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_heat_solver_config.hpp"

namespace lss
{

using memory_space = lss_enumerations::memory_space_enum;
using traverse_direction = lss_enumerations::traverse_direction_enum;
using explicit_pde_schemes = lss_enumerations::explicit_pde_schemes_enum;
using heat_explicit_solver_config_ptr = lss_pde_solvers::heat_explicit_solver_config_ptr;
using heat_explicit_solver_config = lss_pde_solvers::heat_explicit_solver_config;

struct heat_explicit_solver_config_builder
{
  private:
    memory_space memory_space_;
    traverse_direction traverse_direction_;
    explicit_pde_schemes explicit_pde_scheme_;

  public:
    LSS_API explicit heat_explicit_solver_config_builder();

    LSS_API heat_explicit_solver_config_builder &memory(memory_space memory);

    LSS_API heat_explicit_solver_config_builder &traverse(traverse_direction traverse);

    LSS_API heat_explicit_solver_config_builder &explicit_pde_scheme(explicit_pde_schemes explicit_pde_scheme);

    LSS_API heat_explicit_solver_config_ptr build();
};

using namespace lss_pde_solvers::default_heat_solver_configs;

} // namespace lss

#endif ///_HEAT_EXPLICIT_SOLVER_CONFIG_HPP_
