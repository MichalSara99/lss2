#include "lss_grid_config.hpp"

namespace lss_grids
{

grid_config_1d::grid_config_1d(discretization_config_1d_ptr const &discretization_config)
{
    auto const one = static_cast<double>(1.0);
    step_ = one / (discretization_config->number_of_space_points() - 1);
}

double grid_config_1d::step() const
{
    return step_;
}

std::size_t grid_config_1d::index_of(double zeta)
{
    return static_cast<std::size_t>(zeta / step_);
}

double grid_config_1d::value_for(std::size_t idx)
{
    return (step_ * idx);
}

grid_config_2d::grid_config_2d(pde_discretization_config_2d_ptr const &discretization_config)
{
    auto const one = static_cast<double>(1.0);
    step_1_ = one / (discretization_config->number_of_space_points().first - 1);
    step_2_ = one / (discretization_config->number_of_space_points().second - 1);
    grid_1_ = std::make_shared<grid_config_1d>(discretization_config->pde_discretization_1());
    grid_2_ = std::make_shared<grid_config_1d>(discretization_config->pde_discretization_2());
}

grid_config_1d_ptr const &grid_config_2d::grid_1() const
{
    return grid_1_;
};

grid_config_1d_ptr const &grid_config_2d::grid_2() const
{
    return grid_2_;
}

double grid_config_2d::step_1() const
{
    return step_1_;
}

double grid_config_2d::step_2() const
{
    return step_2_;
}

std::size_t grid_config_2d::index_of_1(double zeta)
{
    return static_cast<std::size_t>(zeta / step_1_);
}

std::size_t grid_config_2d::index_of_2(double eta)
{
    return static_cast<std::size_t>(eta / step_2_);
}

double grid_config_2d::value_for_1(std::size_t idx)
{
    return (step_1_ * idx);
}

double grid_config_2d::value_for_2(std::size_t idx)
{
    return (step_2_ * idx);
}

grid_config_3d::grid_config_3d(pde_discretization_config_3d_ptr const &discretization_config)
{
    auto const one = static_cast<double>(1.0);
    step_1_ = one / (std::get<0>(discretization_config->number_of_space_points()) - 1);
    step_2_ = one / (std::get<1>(discretization_config->number_of_space_points()) - 1);
    step_3_ = one / (std::get<2>(discretization_config->number_of_space_points()) - 1);
    grid_1_ = std::make_shared<grid_config_1d>(discretization_config->pde_discretization_1());
    grid_2_ = std::make_shared<grid_config_1d>(discretization_config->pde_discretization_2());
    grid_3_ = std::make_shared<grid_config_1d>(discretization_config->pde_discretization_3());
    grid_12_ = std::make_shared<grid_config_2d>(discretization_config->pde_discretization_12());
    grid_21_ = std::make_shared<grid_config_2d>(discretization_config->pde_discretization_21());
    grid_13_ = std::make_shared<grid_config_2d>(discretization_config->pde_discretization_13());
    grid_31_ = std::make_shared<grid_config_2d>(discretization_config->pde_discretization_31());
    grid_23_ = std::make_shared<grid_config_2d>(discretization_config->pde_discretization_23());
    grid_32_ = std::make_shared<grid_config_2d>(discretization_config->pde_discretization_32());
}

grid_config_1d_ptr const &grid_config_3d::grid_1() const
{
    return grid_1_;
};

grid_config_1d_ptr const &grid_config_3d::grid_2() const
{
    return grid_2_;
}

grid_config_1d_ptr const &grid_config_3d::grid_3() const
{
    return grid_3_;
}

grid_config_2d_ptr const &grid_config_3d::grid_12() const
{
    return grid_12_;
};

grid_config_2d_ptr const &grid_config_3d::grid_21() const
{
    return grid_21_;
};

grid_config_2d_ptr const &grid_config_3d::grid_13() const
{
    return grid_13_;
}

grid_config_2d_ptr const &grid_config_3d::grid_31() const
{
    return grid_31_;
}

grid_config_2d_ptr const &grid_config_3d::grid_23() const
{
    return grid_23_;
}

grid_config_2d_ptr const &grid_config_3d::grid_32() const
{
    return grid_32_;
}

double grid_config_3d::step_1() const
{
    return step_1_;
}

double grid_config_3d::step_2() const
{
    return step_2_;
}

double grid_config_3d::step_3() const
{
    return step_3_;
}

std::size_t grid_config_3d::index_of_1(double zeta)
{
    return static_cast<std::size_t>(zeta / step_1_);
}

std::size_t grid_config_3d::index_of_2(double eta)
{
    return static_cast<std::size_t>(eta / step_2_);
}

std::size_t grid_config_3d::index_of_3(double ny)
{
    return static_cast<std::size_t>(ny / step_3_);
}

double grid_config_3d::value_for_1(std::size_t idx)
{
    return (step_1_ * idx);
}

double grid_config_3d::value_for_2(std::size_t idx)
{
    return (step_2_ * idx);
}

double grid_config_3d::value_for_3(std::size_t idx)
{
    return (step_3_ * idx);
}

} // namespace lss_grids
