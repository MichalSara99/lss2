#include "output.hpp"

namespace lss
{

void output_curve(discretization_config_ptr const &discretization_config, grid_config_1d_ptr const &grid_config_hints,
                  container_t const &container, to to, std::ostream &out)
{
    if (to == to::stream)
    {
        lss_print::print(discretization_config, grid_config_hints, container, out);
    }
    else if (to == to::xml)
    {
        lss_xml::xml(discretization_config, grid_config_hints, container, out);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void output_surface(discretization_config_1d_ptr const &discretization_config,
                    grid_config_1d_ptr const &grid_config_hints, matrix_2d const &container, to to, std::ostream &out)
{
    if (to == to::stream)
    {
        lss_print::print(discretization_config, grid_config_hints, container, out);
    }
    else if (to == to::xml)
    {
        lss_xml::xml(discretization_config, grid_config_hints, container, out);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void output_surface(discretization_config_2d_ptr const &discretization_config,
                    grid_config_2d_ptr const &grid_config_hints, matrix_2d const &container, to to, std::ostream &out)
{
    if (to == to::stream)
    {
        lss_print::print(discretization_config, grid_config_hints, container, out);
    }
    else if (to == to::xml)
    {
        lss_xml::xml(discretization_config, grid_config_hints, container, out);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

void output_surface(discretization_config_2d_ptr const &discretization_config,
                    grid_config_2d_ptr const &grid_config_hints, matrix_3d const &container, to to, std::ostream &out)
{
    if (to == to::stream)
    {
        lss_print::print(discretization_config, grid_config_hints, container, out);
    }
    else if (to == to::xml)
    {
        lss_xml::xml(discretization_config, grid_config_hints, container, out);
    }
    else
    {
        throw std::exception("Unreachable");
    }
}

} // namespace lss
