#if !defined(_ODE_SOLVER_CONFIG_HPP_)
#define _ODE_SOLVER_CONFIG_HPP_

#include "../../../../common/lss_enumerations.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../common/lss_utility.hpp"
#include "../../../../ode_solvers/lss_ode_solver_config.hpp"

namespace lss
{

using memory_space = lss_enumerations::memory_space_enum;
using tridiagonal_method = lss_enumerations::tridiagonal_method_enum;
using factorization = lss_enumerations::factorization_enum;
using ode_implicit_solver_config_ptr = lss_ode_solvers::ode_implicit_solver_config_ptr;
using ode_implicit_solver_config = lss_ode_solvers::ode_implicit_solver_config;

struct ode_implicit_solver_config_builder
{
  private:
    memory_space memory_space_;
    tridiagonal_method tridiagonal_method_;
    factorization tridiagonal_factorization_;

  public:
    LSS_API explicit ode_implicit_solver_config_builder();

    LSS_API ode_implicit_solver_config_builder &memory(memory_space memory);

    LSS_API ode_implicit_solver_config_builder &method(tridiagonal_method method);

    LSS_API ode_implicit_solver_config_builder &tridiagonal_factorization(factorization tridiagonal_factorization);

    LSS_API ode_implicit_solver_config_ptr build();
};

using namespace lss_ode_solvers::default_ode_solver_configs;

} // namespace lss

#endif ///_ODE_SOLVER_CONFIG_HPP_
