#include "lss_xml.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>

//#include "../discretization/lss_grid.hpp"
//#include "../discretization/lss_grid_config.hpp"
//#include "../discretization/lss_grid_transform_config.hpp"

namespace lss_xml
{
// using lss_enumerations::grid_enum;
// using lss_grids::grid_1d;
// using lss_grids::grid_2d;
// using lss_grids::grid_3d;
// using lss_grids::grid_config_1d;
// using lss_grids::grid_config_2d;
// using lss_grids::grid_config_3d;
// using lss_grids::grid_transform_config_1d;
// using lss_grids::grid_transform_config_2d;
// using lss_grids::grid_transform_config_3d;
/*
void xml(discretization_config_1d_ptr const &discretization_config, grid_config_hints_1d_ptr const &grid_hints_cfg,
         container_t const &container, std::ostream &out)
{
    const std::size_t space_size = discretization_config->number_of_space_points();
    LSS_ASSERT(container.size() == space_size, "Container size differs from passed discretization");
    // create grid_config:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(discretization_config);
    const std::string grid_type = (grid_hints_cfg->grid() == grid_enum::Uniform) ? "Uniform" : "Nonuniform";
    // create grid_transform_config:
    auto const &grid_trans_cfg = std::make_shared<grid_transform_config_1d>(discretization_config, grid_hints_cfg);
    double zeta{};

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><CURVE><ATTRIBUTES><GRID>" << grid_type << "</GRID>";
    if (grid_hints_cfg->grid() == grid_enum::Nonuniform)
    {
        out << "<ACCUMULATION_POINT>" << grid_hints_cfg->accumulation_point() << "</ACCUMULATION_POINT>";
    }
    out << "</ATTRIBUTES><AXES><ABSCISSA>";
    out << "<SIZE>" << space_size << "</SIZE><POINTS>";
    for (std::size_t t = 0; t < space_size - 1; ++t)
    {
        zeta = grid_1d::value(grid_cfg, t);
        out << grid_1d::transformed_value(grid_trans_cfg, zeta);
        out << ",";
    }
    zeta = grid_1d::value(grid_cfg, space_size - 1);
    out << grid_1d::transformed_value(grid_trans_cfg, zeta);
    out << "</POINTS></ABSCISSA>";
    out << "<ORDINATE><SIZE>" << space_size << "</SIZE><VALUES>";
    for (std::size_t t = 0; t < space_size - 1; ++t)
    {
        out << container[t];
        out << ",";
    }
    out << container[space_size - 1];
    out << "</VALUES></ORDINATE></AXES></CURVE>";
}

void xml(pde_discretization_config_1d_ptr const &pde_discretization_config,
         grid_config_hints_1d_ptr const &grid_config_hints, container_2d<by_enum::Row> const &container,
         std::ostream &out)
{
    const std::size_t space_size = pde_discretization_config->number_of_space_points();
    const std::size_t time_size = pde_discretization_config->number_of_time_points();
    LSS_ASSERT((container.columns() == space_size) && (container.rows() == time_size),
               "Container size differs from passed discretization");
    // create grid_config:
    auto const &grid_cfg = std::make_shared<grid_config_1d>(pde_discretization_config);
    const std::string grid_type = (grid_config_hints->grid() == grid_enum::Uniform) ? "Uniform" : "Nonuniform";
    // create grid_transform_config:
    auto const &grid_trans_cfg =
        std::make_shared<grid_transform_config_1d>(pde_discretization_config, grid_config_hints);

    const double k = pde_discretization_config->time_step();
    const auto &time = pde_discretization_config->time_range();
    const double time_start = time->lower();
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><SURFACE><ATTRIBUTES><GRID>" << grid_type << "</GRID>";
    out << "<TYPE>SPACE_TIME</TYPE>";
    if (grid_config_hints->grid() == grid_enum::Nonuniform)
    {
        out << "<ACCUMULATION_POINT>" << grid_config_hints->accumulation_point() << "</ACCUMULATION_POINT>";
    }
    out << "</ATTRIBUTES><AXES><ABSCISSA><SIZE>" << space_size << "</SIZE><SPACE_POINTS>";
    double m{}, zeta{};
    for (std::size_t t = 0; t < space_size - 1; ++t)
    {
        zeta = grid_1d::value(grid_cfg, t);
        out << grid_1d::transformed_value(grid_trans_cfg, zeta);
        out << ",";
    }
    zeta = grid_1d::value(grid_cfg, space_size - 1);
    out << grid_1d::transformed_value(grid_trans_cfg, zeta);
    out << "</SPACE_POINTS></ABSCISSA><ABSCISSA><SIZE>" << time_size << "</SIZE><TIME_POINTS>";
    for (std::size_t t = 0; t < time_size - 1; ++t)
    {
        m = static_cast<double>(t);
        out << (time_start + m * k);
        out << ",";
    }
    m = static_cast<double>(time_size - 1);
    out << (time_start + m * k);
    out << "</TIME_POINTS></ABSCISSA><ORDINATE><ROW_SIZE>" << time_size << "</ROW_SIZE><COLUMN_SIZE>" << space_size
        << "</COLUMN_SIZE><VALUES>";
    auto const &values = container.data();
    for (std::size_t t = 0; t < values.size() - 1; ++t)
    {
        out << values[t] << ",";
    }
    out << values[values.size() - 1];
    out << "</VALUES></ORDINATE></AXES></SURFACE>";
}

void xml(pde_discretization_config_2d_ptr const &pde_discretization_config,
         grid_config_hints_2d_ptr const &grid_config_hints, container_2d<by_enum::Row> const &container,
         std::ostream &out)
{
    const auto &space_sizes = pde_discretization_config->number_of_space_points();
    LSS_ASSERT((container.columns() == std::get<1>(space_sizes)) && (container.rows() == std::get<0>(space_sizes)),
               "The input cont container must have the correct size");
    // create grid_config:
    auto const &grid_cfg = std::make_shared<grid_config_2d>(pde_discretization_config);
    // create grid_transform_config:
    auto const &grid_trans_cfg =
        std::make_shared<grid_transform_config_2d>(pde_discretization_config, grid_config_hints);
    const std::string grid_type = (grid_config_hints->grid() == grid_enum::Uniform) ? "Uniform" : "Nonuniform";
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><SURFACE><ATTRIBUTES><GRID>" << grid_type << "</GRID>";
    out << "<TYPE>SPACE_SPACE</TYPE>";
    if (grid_config_hints->grid() == grid_enum::Nonuniform)
    {
        out << "<ACCUMULATION_POINT>" << grid_config_hints->accumulation_point() << "</ACCUMULATION_POINT>";
    }
    out << "</ATTRIBUTES><AXES><ABSCISSA><SIZE>" << space_sizes.first << "</SIZE><POINTS>";

    double zeta{}, eta{};
    for (std::size_t t = 0; t < space_sizes.first - 1; ++t)
    {
        zeta = grid_2d::value_1(grid_cfg, t);
        out << grid_2d::transformed_value_1(grid_trans_cfg, zeta);
        out << ",";
    }
    zeta = grid_2d::value_1(grid_cfg, space_sizes.first - 1);
    out << grid_2d::transformed_value_1(grid_trans_cfg, zeta);
    out << "</POINTS></ABSCISSA><ABSCISSA><SIZE>" << space_sizes.second << "</SIZE><POINTS>";
    for (std::size_t t = 0; t < space_sizes.second - 1; ++t)
    {
        eta = grid_2d::value_2(grid_cfg, t);
        out << grid_2d::transformed_value_2(grid_trans_cfg, eta);
        out << ",";
    }
    eta = grid_2d::value_2(grid_cfg, space_sizes.second - 1);
    out << grid_2d::transformed_value_2(grid_trans_cfg, eta);
    out << "</POINTS></ABSCISSA><ORDINATE><ROW_SIZE>" << space_sizes.first << "</ROW_SIZE><COLUMN_SIZE>"
        << space_sizes.second << "</COLUMN_SIZE><VALUES>";
    auto const &values = container.data();
    for (std::size_t t = 0; t < values.size() - 1; ++t)
    {
        out << values[t] << ",";
    }
    out << values[values.size() - 1];
    out << "</VALUES></ORDINATE></AXES></SURFACE>";
}

void xml(pde_discretization_config_3d_ptr const &pde_discretization_config,
         grid_config_hints_3d_ptr const &grid_config_hints, container_3d<by_enum::LayerPlane> const &container,
         std::ostream &out)
{
    const auto &space_sizes = pde_discretization_config->number_of_space_points();
    LSS_ASSERT((container.columns() == std::get<1>(space_sizes)) && (container.rows() == std::get<0>(space_sizes)) &&
                   (container.layers() == std::get<2>(space_sizes)),
               "The input cont container must have the correct size");
    // create grid_config:
    auto const &grid_cfg = std::make_shared<grid_config_3d>(pde_discretization_config);
    // create grid_transform_config:
    auto const &grid_trans_cfg =
        std::make_shared<grid_transform_config_3d>(pde_discretization_config, grid_config_hints);
    const std::string grid_type = (grid_config_hints->grid() == grid_enum::Uniform) ? "Uniform" : "Nonuniform";
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><SURFACE><ATTRIBUTES><GRID>" << grid_type << "</GRID>";
    out << "<TYPE>SPACE_SPACE_SPACE</TYPE>";
    if (grid_config_hints->grid() == grid_enum::Nonuniform)
    {
        out << "<ACCUMULATION_POINT>" << grid_config_hints->accumulation_point() << "</ACCUMULATION_POINT>";
    }
    const std::size_t space_size_x = std::get<0>(space_sizes);
    out << "</ATTRIBUTES><AXES><ABSCISSA><SIZE>" << space_size_x << "</SIZE><POINTS>";
    double zeta{}, eta{}, ny{};
    for (std::size_t t = 0; t < space_size_x - 1; ++t)
    {
        zeta = grid_3d::value_1(grid_cfg, t);
        out << grid_3d::transformed_value_1(grid_trans_cfg, zeta);
        out << ",";
    }
    zeta = grid_3d::value_1(grid_cfg, space_size_x - 1);
    out << grid_3d::transformed_value_1(grid_trans_cfg, zeta);
    const std::size_t space_size_y = std::get<1>(space_sizes);
    out << "</POINTS></ABSCISSA><ABSCISSA><SIZE>" << space_size_y << "</SIZE><POINTS>";
    for (std::size_t t = 0; t < space_size_y - 1; ++t)
    {
        eta = grid_3d::value_2(grid_cfg, t);
        out << grid_3d::transformed_value_2(grid_trans_cfg, eta);
        out << ",";
    }
    eta = grid_3d::value_2(grid_cfg, space_size_y - 1);
    out << grid_3d::transformed_value_2(grid_trans_cfg, eta);
    const std::size_t space_size_z = std::get<2>(space_sizes);
    out << "</POINTS></ABSCISSA><ABSCISSA><SIZE>" << space_size_z << "</SIZE><POINTS>";
    for (std::size_t t = 0; t < space_size_z - 1; ++t)
    {
        ny = grid_3d::value_3(grid_cfg, t);
        out << grid_3d::transformed_value_3(grid_trans_cfg, ny);
        out << ",";
    }
    ny = grid_3d::value_3(grid_cfg, space_size_z - 1);
    out << grid_3d::transformed_value_3(grid_trans_cfg, ny);

    out << "</POINTS></ABSCISSA><ORDINATE><ROW_SIZE>" << space_size_x << "</ROW_SIZE><COLUMN_SIZE>" << space_size_y
        << "</COLUMN_SIZE><LAYER_SIZE>" << space_size_z << "</LAYER_SIZE><VALUES>";
    auto const &values = container.data();
    for (std::size_t t = 0; t < values.size() - 1; ++t)
    {
        out << values[t] << ",";
    }
    out << values[values.size() - 1];
    out << "</VALUES></ORDINATE></AXES></SURFACE>";
}
*/
} // namespace lss_xml
