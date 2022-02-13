#include "explicit_heat_solver_1d.hpp"

namespace lss
{

explicit_heat_solver_1d_builder::explicit_heat_solver_1d_builder()
{
}

explicit_heat_solver_1d_builder &explicit_heat_solver_1d_builder::heat_data_config(
    const heat_data_config_1d_ptr &heat_data_config)
{
    heat_data_config_ = heat_data_config;
    return *this;
}

explicit_heat_solver_1d_builder &explicit_heat_solver_1d_builder::discretization_config(
    const discretization_config_1d_ptr &discretization_config)
{
    discretization_config_ = discretization_config;
    return *this;
}

explicit_heat_solver_1d_builder &explicit_heat_solver_1d_builder::boundary_pair(const boundary_1d_pair &boundary_pair)
{
    boundary_pair_ = boundary_pair;
    return *this;
}

explicit_heat_solver_1d_builder &explicit_heat_solver_1d_builder::grid_hints(const grid_config_1d_ptr &grid_hints)
{
    grid_cfg_ = grid_hints;
    return *this;
}

explicit_heat_solver_1d_builder &explicit_heat_solver_1d_builder::solver_config(
    const heat_explicit_solver_config_ptr &solver_config)
{
    solver_config_ = solver_config;
    return *this;
}

explicit_heat_solver_1d_ptr explicit_heat_solver_1d_builder::build()
{
    return std::make_shared<explicit_heat_solver_1d>(heat_data_config_, discretization_config_, boundary_pair_,
                                                     grid_cfg_, solver_config_);
}

} // namespace lss
