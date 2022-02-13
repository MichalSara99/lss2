#include "wave_implicit_solver_config.hpp"

namespace lss
{

wave_implicit_solver_config_builder::wave_implicit_solver_config_builder()
{
}

wave_implicit_solver_config_builder &wave_implicit_solver_config_builder::memory(memory_space memory)
{
    memory_space_ = memory;
    return *this;
}

wave_implicit_solver_config_builder &wave_implicit_solver_config_builder::traverse(traverse_direction traverse)
{
    traverse_direction_ = traverse;
    return *this;
}

wave_implicit_solver_config_builder &wave_implicit_solver_config_builder::method(tridiagonal_method method)
{
    tridiagonal_method_ = method;
    return *this;
}

wave_implicit_solver_config_builder &wave_implicit_solver_config_builder::tridiagonal_factorization(
    factorization tridiagonal_factorization)
{
    tridiagonal_factorization_ = tridiagonal_factorization;
    return *this;
}

wave_implicit_solver_config_ptr wave_implicit_solver_config_builder::build()
{
    return std::make_shared<wave_implicit_solver_config>(memory_space_, traverse_direction_, tridiagonal_method_,
                                                         tridiagonal_factorization_);
}
} // namespace lss
