#include "lss_grid_transform_config.hpp"

namespace lss_grids
{

grid_transform_config_1d::grid_transform_config_1d(discretization_config_1d_ptr const &discretization_config,
                                                   grid_config_hints_1d_ptr const &grid_hints)
{
    auto const low = discretization_config->space_range()->lower();
    auto const high = discretization_config->space_range()->upper();
    initialize(low, high, grid_hints);
}

void grid_transform_config_1d::initialize(double low, double high, grid_config_hints_1d_ptr const &grid_hints)
{
    auto const point = grid_hints->accumulation_point();
    alpha_ = high - low;
    // in case non-uniform spacing is requested alpha and beta are overriden
    if (grid_hints->grid() == grid_enum::Nonuniform)
    {
        alpha_ = (point - low) / grid_hints->alpha_scale();
    }
    init_ = point;
    c_[0] = std::asinh((high - point) / alpha_);
    c_[1] = std::asinh((low - point) / alpha_);
    auto const c_diff = (c_[0] - c_[1]);

    // initialize derivatives:
    auto const one = static_cast<double>(1.0);
    a_der_ = [=](double zeta) { return (alpha_ * c_diff * std::cosh(c_[0] * zeta + c_[1] * (one - zeta))); };
    b_der_ = [=](double zeta) { return (alpha_ * c_diff * c_diff * std::sinh(c_[0] * zeta + c_[1] * (one - zeta))); };
}

std::function<double(double)> const &grid_transform_config_1d::a_derivative() const
{
    return a_der_;
}

std::function<double(double)> const &grid_transform_config_1d::b_derivative() const
{
    return b_der_;
}

double grid_transform_config_1d::value_for(double zeta)
{
    auto const one = static_cast<double>(1.0);
    return (init_ + alpha_ * std::sinh(c_[0] * zeta + c_[1] * (one - zeta)));
}

grid_transform_config_2d::grid_transform_config_2d(pde_discretization_config_2d_ptr const &discretization_config,
                                                   grid_config_hints_2d_ptr const &grid_hints)
{
    initialize(discretization_config, grid_hints);
}

void grid_transform_config_2d::initialize(pde_discretization_config_2d_ptr const &discretization_config,
                                          grid_config_hints_2d_ptr const &grid_hints)
{
    auto const low_1 = discretization_config->space_range().first->lower();
    auto const high_1 = discretization_config->space_range().first->upper();
    auto const low_2 = discretization_config->space_range().second->lower();
    auto const high_2 = discretization_config->space_range().second->upper();
    auto const point = grid_hints->accumulation_point();
    alpha_ = high_1 - low_1;
    beta_ = high_2 - low_2;
    // in case non-uniform spacing is requested alpha and beta are overriden
    if (grid_hints->grid() == grid_enum::Nonuniform)
    {
        alpha_ = (point - low_1) / grid_hints->alpha_scale();
        beta_ = (high_2 - low_2) / grid_hints->beta_scale();
    }
    init_1_ = point;
    init_2_ = low_2;
    c_[0] = std::asinh((high_1 - point) / alpha_);
    c_[1] = std::asinh((low_1 - point) / alpha_);
    d_ = std::asinh((high_2 - low_2) / beta_);
    auto const c_diff = (c_[0] - c_[1]);

    // initialize derivatives:
    auto const one = static_cast<double>(1.0);
    a_der_ = [=](double zeta) { return (alpha_ * c_diff * std::cosh(c_[0] * zeta + c_[1] * (one - zeta))); };
    b_der_ = [=](double eta) { return (beta_ * d_ * std::cosh(d_ * eta)); };
    c_der_ = [=](double zeta) { return (alpha_ * c_diff * c_diff * std::sinh(c_[0] * zeta + c_[1] * (one - zeta))); };
    d_der_ = [=](double eta) { return (beta_ * d_ * d_ * std::sinh(d_ * eta)); };
}

std::function<double(double)> const &grid_transform_config_2d::a_derivative() const
{
    return a_der_;
}

std::function<double(double)> const &grid_transform_config_2d::b_derivative() const
{
    return b_der_;
}

std::function<double(double)> const &grid_transform_config_2d::c_derivative() const
{
    return c_der_;
}

std::function<double(double)> const &grid_transform_config_2d::d_derivative() const
{
    return d_der_;
}

double grid_transform_config_2d::value_for_1(double zeta)
{
    auto const one = static_cast<double>(1.0);
    return (init_1_ + alpha_ * std::sinh(c_[0] * zeta + c_[1] * (one - zeta)));
}

double grid_transform_config_2d::value_for_2(double eta)
{
    return (init_2_ + beta_ * std::sinh(d_ * eta));
}

grid_transform_config_3d::grid_transform_config_3d(pde_discretization_config_3d_ptr const &discretization_config,
                                                   grid_config_hints_3d_ptr const &grid_hints)
{
    initialize(discretization_config, grid_hints);
}

void grid_transform_config_3d::initialize(pde_discretization_config_3d_ptr const &discretization_config,
                                          grid_config_hints_3d_ptr const &grid_hints)
{

    auto const low_1 = std::get<0>(discretization_config->space_range())->lower();
    auto const high_1 = std::get<0>(discretization_config->space_range())->upper();
    auto const low_2 = std::get<1>(discretization_config->space_range())->lower();
    auto const high_2 = std::get<1>(discretization_config->space_range())->upper();
    auto const low_3 = std::get<2>(discretization_config->space_range())->lower();
    auto const high_3 = std::get<2>(discretization_config->space_range())->upper();
    auto const point = grid_hints->accumulation_point();
    alpha_ = high_1 - low_1;
    beta_1_ = high_2 - low_2;
    beta_2_ = high_3 - low_3;
    // in case non-uniform spacing is requested alpha and beta are overriden
    if (grid_hints->grid() == grid_enum::Nonuniform)
    {
        alpha_ = (point - low_1) / grid_hints->alpha_scale();
        beta_1_ = (high_2 - low_2) / std::get<0>(grid_hints->beta_scales());
        beta_2_ = (high_3 - low_3) / std::get<0>(grid_hints->beta_scales());
    }
    init_1_ = point;
    init_2_ = low_2;
    init_3_ = low_3;
    c_[0] = std::asinh((high_1 - point) / alpha_);
    c_[1] = std::asinh((low_1 - point) / alpha_);
    d_ = std::asinh((high_2 - low_2) / beta_1_);
    e_ = std::asinh((high_3 - low_3) / beta_2_);
    auto const c_diff = (c_[0] - c_[1]);

    // initialize derivatives:
    auto const one = static_cast<double>(1.0);
    a_1_der_ = [=](double zeta) { return (alpha_ * c_diff * std::cosh(c_[0] * zeta + c_[1] * (one - zeta))); };
    a_2_der_ = [=](double zeta) { return (alpha_ * c_diff * c_diff * std::sinh(c_[0] * zeta + c_[1] * (one - zeta))); };

    b_1_der_ = [=](double eta) { return (beta_1_ * d_ * std::cosh(d_ * eta)); };
    b_2_der_ = [=](double eta) { return (beta_1_ * d_ * d_ * std::sinh(d_ * eta)); };

    c_1_der_ = [=](double ny) { return (beta_2_ * e_ * std::cosh(e_ * ny)); };
    c_2_der_ = [=](double ny) { return (beta_2_ * e_ * e_ * std::sinh(e_ * ny)); };
}

std::function<double(double)> const &grid_transform_config_3d::a_1_derivative() const
{
    return a_1_der_;
}

std::function<double(double)> const &grid_transform_config_3d::a_2_derivative() const
{
    return a_2_der_;
}

std::function<double(double)> const &grid_transform_config_3d::b_1_derivative() const
{
    return b_1_der_;
}

std::function<double(double)> const &grid_transform_config_3d::b_2_derivative() const
{
    return b_2_der_;
}

std::function<double(double)> const &grid_transform_config_3d::c_1_derivative() const
{
    return c_1_der_;
}

std::function<double(double)> const &grid_transform_config_3d::c_2_derivative() const
{
    return c_2_der_;
}

double grid_transform_config_3d::value_for_1(double zeta)
{
    auto const one = static_cast<double>(1.0);
    return (init_1_ + alpha_ * std::sinh(c_[0] * zeta + c_[1] * (one - zeta)));
}

double grid_transform_config_3d::value_for_2(double eta)
{
    return (init_2_ + beta_1_ * std::sinh(d_ * eta));
}

double grid_transform_config_3d::value_for_3(double ny)
{
    return (init_3_ + beta_2_ * std::sinh(e_ * ny));
}

} // namespace lss_grids
