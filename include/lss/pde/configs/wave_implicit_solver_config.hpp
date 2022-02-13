#if !defined(_WAVE_IMPLICIT_SOLVER_CONFIG_HPP_)
#define _WAVE_IMPLICIT_SOLVER_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../pde_solvers/lss_wave_solver_config.hpp"

namespace lss
{

using memory_space = lss_enumerations::memory_space_enum;
using traverse_direction = lss_enumerations::traverse_direction_enum;
using tridiagonal_method = lss_enumerations::tridiagonal_method_enum;
using factorization = lss_enumerations::factorization_enum;
using wave_implicit_solver_config_ptr = lss_pde_solvers::wave_implicit_solver_config_ptr;
using wave_implicit_solver_config = lss_pde_solvers::wave_implicit_solver_config;

struct wave_implicit_solver_config_builder
{
  private:
    memory_space memory_space_;
    traverse_direction traverse_direction_;
    tridiagonal_method tridiagonal_method_;
    factorization tridiagonal_factorization_;

  public:
    LSS_API explicit wave_implicit_solver_config_builder();

    LSS_API wave_implicit_solver_config_builder &memory(memory_space memory);

    LSS_API wave_implicit_solver_config_builder &traverse(traverse_direction traverse);

    LSS_API wave_implicit_solver_config_builder &method(tridiagonal_method method);

    LSS_API wave_implicit_solver_config_builder &tridiagonal_factorization(factorization tridiagonal_factorization);

    LSS_API wave_implicit_solver_config_ptr build();
};

using namespace lss_pde_solvers::default_wave_solver_configs;

} // namespace lss

#endif ///_WAVE_IMPLICIT_SOLVER_CONFIG_HPP_
