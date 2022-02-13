#include "implicit_wave_solver_1d.hpp"

namespace lss
{

implicit_wave_solver_1d_builder::implicit_wave_solver_1d_builder()
{
}

implicit_wave_solver_1d_builder &implicit_wave_solver_1d_builder::heat_data_config(
    const wave_data_config_1d_ptr &wave_data_config)
{
    wave_data_config_ = wave_data_config;
    return *this;
}

implicit_wave_solver_1d_builder &implicit_wave_solver_1d_builder::discretization_config(
    const discretization_config_1d_ptr &discretization_config)
{
    discretization_config_ = discretization_config;
    return *this;
}

implicit_wave_solver_1d_builder &implicit_wave_solver_1d_builder::boundary_pair(const boundary_1d_pair &boundary_pair)
{
    boundary_pair_ = boundary_pair;
    return *this;
}

implicit_wave_solver_1d_builder &implicit_wave_solver_1d_builder::grid_hints(const grid_config_1d_ptr &grid_hints)
{
    grid_cfg_ = grid_hints;
    return *this;
}

implicit_wave_solver_1d_builder &implicit_wave_solver_1d_builder::solver_config(
    const wave_implicit_solver_config_ptr &solver_config)
{
    solver_config_ = solver_config;
    return *this;
}

implicit_wave_solver_1d_builder &implicit_wave_solver_1d_builder::solver_config_details(
    const std::map<std::string, double> &solver_config_details)
{
    solver_config_details_ = solver_config_details;
    return *this;
}

implicit_wave_solver_1d_ptr implicit_wave_solver_1d_builder::build()
{
    return std::make_shared<implicit_wave_solver_1d>(wave_data_config_, discretization_config_, boundary_pair_,
                                                     grid_cfg_, solver_config_, solver_config_details_);
}

} // namespace lss
