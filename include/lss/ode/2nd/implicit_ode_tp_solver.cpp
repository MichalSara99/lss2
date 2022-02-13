#include "implicit_ode_tp_solver.hpp"

namespace lss
{

implicit_ode_tp_solver_builder::implicit_ode_tp_solver_builder()
{
}

implicit_ode_tp_solver_builder &implicit_ode_tp_solver_builder::ode_data_config(
    const ode_data_config_ptr &ode_data_config)
{
    ode_data_cfg_ = ode_data_config;
    return *this;
}

implicit_ode_tp_solver_builder &implicit_ode_tp_solver_builder::discretization_config(
    const discretization_config_ode_ptr &discretization_config)
{
    discretization_cfg_ = discretization_config;
    return *this;
}

implicit_ode_tp_solver_builder &implicit_ode_tp_solver_builder::boundary_pair(const boundary_1d_pair &boundary_pair)
{
    boundary_pair_ = boundary_pair;
    return *this;
}

implicit_ode_tp_solver_builder &implicit_ode_tp_solver_builder::grid_hints(const grid_config_1d_ptr &grid_hints)
{
    grid_cfg_ = grid_hints;
    return *this;
}

implicit_ode_tp_solver_builder &implicit_ode_tp_solver_builder::solver_config(
    const ode_implicit_solver_config_ptr &solver_config)
{
    ode_solver_config_ = solver_config;
    return *this;
}

implicit_ode_tp_solver_builder &implicit_ode_tp_solver_builder::solver_config_details(
    const std::map<std::string, double> &solver_config_details)
{
    ode_solver_config_details_ = solver_config_details;
    return *this;
}

implicit_ode_tp_solver_ptr implicit_ode_tp_solver_builder::build()
{
    return std::make_shared<implicit_ode_tp_solver>(ode_data_cfg_, discretization_cfg_, boundary_pair_, grid_cfg_,
                                                    ode_solver_config_, ode_solver_config_details_);
}

} // namespace lss
