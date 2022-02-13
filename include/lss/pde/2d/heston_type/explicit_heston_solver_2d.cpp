#include "explicit_heston_solver_2d.hpp"

namespace lss
{

explicit_heston_solver_2d_builder::explicit_heston_solver_2d_builder()
{
}

explicit_heston_solver_2d_builder &explicit_heston_solver_2d_builder::heat_data_config(
    const heat_data_config_2d_ptr &heat_data_config)
{
    heat_data_config_ = heat_data_config;
    return *this;
}

explicit_heston_solver_2d_builder &explicit_heston_solver_2d_builder::discretization_config(
    const discretization_config_2d_ptr &discretization_config)
{
    discretization_config_ = discretization_config;
    return *this;
}

explicit_heston_solver_2d_builder &explicit_heston_solver_2d_builder::vertical_upper_boundary(
    const boundary_2d_ptr &boundary_ptr)
{
    vertical_upper_boundary_ = boundary_ptr;
    return *this;
}

explicit_heston_solver_2d_builder &explicit_heston_solver_2d_builder::horizontal_boundary_pair(
    const boundary_2d_pair &boundary_pair)
{
    horizontal_boundary_pair_ = boundary_pair;
    return *this;
}

explicit_heston_solver_2d_builder &explicit_heston_solver_2d_builder::grid_hints(const grid_config_2d_ptr &grid_hints)
{
    grid_cfg_ = grid_hints;
    return *this;
}

explicit_heston_solver_2d_builder &explicit_heston_solver_2d_builder::solver_config(
    const heat_explicit_solver_config_ptr &solver_config)
{
    solver_config_ = solver_config;
    return *this;
}

explicit_heston_solver_2d_ptr explicit_heston_solver_2d_builder::build()
{
    return std::make_shared<explicit_heston_solver_2d>(heat_data_config_, discretization_config_,
                                                       vertical_upper_boundary_, horizontal_boundary_pair_, grid_cfg_,
                                                       solver_config_);
}

} // namespace lss
