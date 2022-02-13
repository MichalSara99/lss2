#include "lss_discretization.hpp"
#include "lss_grid.hpp"

using lss_grids::grid_1d;
using lss_grids::grid_2d;
using lss_grids::grid_3d;
using lss_grids::grid_config_1d_ptr;
using lss_grids::grid_config_2d_ptr;
using lss_grids::grid_config_3d_ptr;

// discretization_1d

void discretization_1d::of_space(grid_config_1d_ptr const &grid_config, container_t &output)
{
    LSS_ASSERT(output.size() > 0, "The output container must be initialized.");
    for (std::size_t t = 0; t < output.size(); ++t)
    {
        output[t] = grid_1d::value(grid_config, t);
    }
}

void discretization_1d::of_function(grid_config_1d_ptr const &grid_config, double const &time,
                                    std::function<double(double, double)> const &fun, container_t &output)
{
    LSS_ASSERT(output.size() > 0, "The output container must be initialized.");
    for (std::size_t t = 0; t < output.size(); ++t)
    {
        output[t] = fun(time, grid_1d::value(grid_config, t));
    }
}

void discretization_1d::of_function(grid_config_1d_ptr const &grid_config, std::function<double(double)> const &fun,
                                    container_t &output)
{
    LSS_ASSERT(output.size() > 0, "The output container must be initialized.");
    for (std::size_t t = 0; t < output.size(); ++t)
    {
        output[t] = fun(grid_1d::value(grid_config, t));
    }
}

// discretization_2d

void discretization_2d::of_function(grid_config_2d_ptr const &grid_config,
                                    std::function<double(double, double)> const &fun, matrix_2d &output)
{
    LSS_ASSERT(output.rows() > 0, "The output container must be initialized.");
    LSS_ASSERT(output.columns() > 0, "The output container must be initialized.");
    for (std::size_t r = 0; r < output.rows(); ++r)
    {
        for (std::size_t c = 0; c < output.columns(); ++c)
        {
            output(r, c) = fun(grid_2d::value_1(grid_config, r), grid_2d::value_2(grid_config, c));
        }
    }
}

void discretization_2d::of_function(grid_config_2d_ptr const &grid_config, double const &time,
                                    std::function<double(double, double, double)> const &fun, matrix_2d &output)
{
    LSS_ASSERT(output.rows() > 0, "The output container must be initialized.");
    LSS_ASSERT(output.columns() > 0, "The output container must be initialized.");
    for (std::size_t r = 0; r < output.rows(); ++r)
    {
        for (std::size_t c = 0; c < output.columns(); ++c)
        {
            output(r, c) = fun(time, grid_2d::value_1(grid_config, r), grid_2d::value_2(grid_config, c));
        }
    }
}

void discretization_2d::of_function(grid_config_2d_ptr const &grid_config, double const &time,
                                    std::function<double(double, double, double)> const &fun, std::size_t const &rows,
                                    std::size_t const &cols, container_t &cont)
{
    LSS_ASSERT(cont.size() > 0, "The output container must be initialized.");
    for (std::size_t r = 0; r < rows; ++r)
    {
        for (std::size_t c = 0; c < cols; ++c)
        {
            cont[c + r * cols] = fun(time, grid_2d::value_1(grid_config, r), grid_2d::value_2(grid_config, c));
        }
    }
}

// discretization_3d

void discretization_3d::of_function(grid_config_3d_ptr const &grid_config,
                                    std::function<double(double, double, double)> const &fun, matrix_3d &output)
{
    LSS_ASSERT(output.rows() > 0, "The output container must be initialized.");
    LSS_ASSERT(output.columns() > 0, "The output container must be initialized.");
    LSS_ASSERT(output.layers() > 0, "The output container must be initialized.");
    for (std::size_t l = 0; l < output.layers(); ++l)
    {
        for (std::size_t r = 0; r < output.rows(); ++r)
        {
            for (std::size_t c = 0; c < output.columns(); ++c)
            {
                output(r, c, l) = fun(grid_3d::value_1(grid_config, r), grid_3d::value_2(grid_config, c),
                                      grid_3d::value_3(grid_config, l));
            }
        }
    }
}

void discretization_3d::of_function(grid_config_3d_ptr const &grid_config, double const &time,
                                    std::function<double(double, double, double, double)> const &fun, matrix_3d &output)
{
    LSS_ASSERT(output.rows() > 0, "The output container must be initialized.");
    LSS_ASSERT(output.columns() > 0, "The output container must be initialized.");
    LSS_ASSERT(output.layers() > 0, "The output container must be initialized.");

    for (std::size_t l = 0; l < output.layers(); ++l)
    {
        for (std::size_t r = 0; r < output.rows(); ++r)
        {
            for (std::size_t c = 0; c < output.columns(); ++c)
            {
                output(r, c, l) = fun(time, grid_3d::value_1(grid_config, r), grid_3d::value_2(grid_config, c),
                                      grid_3d::value_3(grid_config, l));
            }
        }
    }
}
