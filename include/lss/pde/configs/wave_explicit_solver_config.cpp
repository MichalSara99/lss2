#include "wave_explicit_solver_config.hpp"

namespace lss
{

wave_explicit_solver_config_builder::wave_explicit_solver_config_builder()
{
}

wave_explicit_solver_config_builder &wave_explicit_solver_config_builder::memory(memory_space memory)
{
    memory_space_ = memory;
    return *this;
}

wave_explicit_solver_config_builder &wave_explicit_solver_config_builder::traverse(traverse_direction traverse)
{
    traverse_direction_ = traverse;
    return *this;
}

wave_explicit_solver_config_ptr wave_explicit_solver_config_builder::build()
{
    return std::make_shared<wave_explicit_solver_config>(memory_space_, traverse_direction_);
}
} // namespace lss
