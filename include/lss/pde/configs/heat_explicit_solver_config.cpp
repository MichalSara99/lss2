#include "heat_explicit_solver_config.hpp"

namespace lss
{

heat_explicit_solver_config_builder::heat_explicit_solver_config_builder()
{
}

heat_explicit_solver_config_builder &heat_explicit_solver_config_builder::memory(memory_space memory)
{
    memory_space_ = memory;
    return *this;
}

heat_explicit_solver_config_builder &heat_explicit_solver_config_builder::traverse(traverse_direction traverse)
{
    traverse_direction_ = traverse;
    return *this;
}

heat_explicit_solver_config_builder &heat_explicit_solver_config_builder::explicit_pde_scheme(
    explicit_pde_schemes explicit_pde_scheme)
{
    explicit_pde_scheme_ = explicit_pde_scheme;
    return *this;
}

heat_explicit_solver_config_ptr heat_explicit_solver_config_builder::build()
{
    return std::make_shared<heat_explicit_solver_config>(memory_space_, traverse_direction_, explicit_pde_scheme_);
}
} // namespace lss
