#include "lss_hhw_boundary_transform.hpp"

#include "../../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../../common/lss_macros.hpp"
#include "../../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_3d;
using lss_boundary::neumann_boundary_3d;
using lss_grids::grid_1d;
using lss_grids::grid_2d;
using lss_grids::grid_3d;

void hhw_boundary_transform::initialize(boundary_3d_pair const &x_boundary_pair,
                                        boundary_3d_ptr const &y_upper_boundary_ptr,
                                        boundary_3d_pair const &z_boundary_pair,
                                        grid_transform_config_3d_ptr const grid_transform_config)
{
    LSS_VERIFY(grid_transform_config, "grid_transform_config must not be null");
    auto const one = 1.0;
    auto const &y_upper_orig = y_upper_boundary_ptr;
    auto const &x_lower_orig = std::get<0>(x_boundary_pair);
    auto const &x_upper_orig = std::get<1>(x_boundary_pair);
    auto const &z_lower_orig = std::get<0>(z_boundary_pair);
    auto const &z_upper_orig = std::get<1>(z_boundary_pair);
    auto const &a_prime = grid_transform_config->a_1_derivative();
    auto const &c_prime = grid_transform_config->c_1_derivative();

    // transform Y's upper:
    auto const y_upper_trans = [=](double t, double zeta, double ny) -> double {
        auto const x = grid_3d::transformed_value_1(grid_transform_config, zeta);
        auto const z = grid_3d::transformed_value_3(grid_transform_config, ny);
        return y_upper_orig->value(t, x, z);
    };
    // transform both Xs:
    // X lower:
    auto const x_lower_trans = [=](double t, double eta, double ny) -> double {
        auto const y = grid_3d::transformed_value_2(grid_transform_config, eta);
        auto const z = grid_3d::transformed_value_3(grid_transform_config, ny);
        return x_lower_orig->value(t, y, z);
    };
    // X upper:
    auto const x_upper_trans = [=](double t, double eta, double ny) -> double {
        auto const y = grid_3d::transformed_value_2(grid_transform_config, eta);
        auto const z = grid_3d::transformed_value_3(grid_transform_config, ny);
        return (x_upper_orig->value(t, y, z) * a_prime(one));
    };

    // transform both Zs:
    // Z lower:
    auto const z_lower_trans = [=](double t, double zeta, double eta) -> double {
        auto const x = grid_3d::transformed_value_1(grid_transform_config, zeta);
        auto const y = grid_3d::transformed_value_2(grid_transform_config, eta);
        return (z_lower_orig->value(t, x, y) * c_prime(0.0));
    };
    // Z upper:
    auto const z_upper_trans = [=](double t, double zeta, double eta) -> double {
        auto const x = grid_3d::transformed_value_1(grid_transform_config, zeta);
        auto const y = grid_3d::transformed_value_2(grid_transform_config, eta);
        return (z_upper_orig->value(t, x, y) * c_prime(one));
    };
    y_upper_ptr_ = std::make_shared<dirichlet_boundary_3d>(y_upper_trans);
    x_pair_ptr_ = std::make_pair(std::make_shared<dirichlet_boundary_3d>(x_lower_trans),
                                 std::make_shared<neumann_boundary_3d>(x_upper_trans));
    z_pair_ptr_ = std::make_pair(std::make_shared<neumann_boundary_3d>(z_lower_trans),
                                 std::make_shared<neumann_boundary_3d>(z_upper_trans));
}

hhw_boundary_transform::hhw_boundary_transform(boundary_3d_pair const &x_boundary_pair,
                                               boundary_3d_ptr const &y_upper_boundary_ptr,
                                               boundary_3d_pair const &z_boundary_pair,
                                               grid_transform_config_3d_ptr const grid_transform_config)
{
    initialize(x_boundary_pair, y_upper_boundary_ptr, z_boundary_pair, grid_transform_config);
}

hhw_boundary_transform::~hhw_boundary_transform()
{
}

boundary_3d_ptr const &hhw_boundary_transform::y_upper_boundary() const
{
    return y_upper_ptr_;
}

boundary_3d_pair const &hhw_boundary_transform::x_boundary_pair() const
{
    return x_pair_ptr_;
}

boundary_3d_pair const &hhw_boundary_transform::z_boundary_pair() const
{
    return z_pair_ptr_;
}

} // namespace lss_pde_solvers
