#include "lss_heston_boundary_transform.hpp"

#include "../../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../../common/lss_macros.hpp"
#include "../../../../discretization/lss_grid.hpp"

namespace lss_pde_solvers
{

using lss_boundary::dirichlet_boundary_2d;
using lss_boundary::neumann_boundary_2d;
using lss_grids::grid_1d;
using lss_grids::grid_2d;

void heston_boundary_transform::initialize(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                           boundary_2d_pair const &horizontal_boundary_pair,
                                           grid_transform_config_2d_ptr const grid_transform_config)
{
    LSS_VERIFY(grid_transform_config, "grid_transform_config must not be null");
    auto const one = 1.0;
    auto const &v_upper_orig = vertical_upper_boundary_ptr;
    auto const &h_lower_orig = std::get<0>(horizontal_boundary_pair);
    auto const &h_upper_orig = std::get<1>(horizontal_boundary_pair);
    auto const &a_der = grid_transform_config->a_derivative();
    std::function<double(double, double)> v_upper_trans;
    std::function<double(double, double)> h_upper_trans;

    if (auto ver_ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(v_upper_orig))
    {
        // transform vertical upper:
        v_upper_trans = [=](double t, double zeta) -> double {
            auto const x = grid_2d::transformed_value_1(grid_transform_config, zeta);
            return v_upper_orig->value(t, x);
        };
        v_upper_ptr_ = std::make_shared<dirichlet_boundary_2d>(v_upper_trans);
    }
    else if (auto ver_ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(v_upper_orig))
    {
        // transform vertical upper:
        auto const &b_der = grid_transform_config->b_derivative();
        v_upper_trans = [=](double t, double zeta) -> double {
            auto const x = grid_2d::transformed_value_1(grid_transform_config, zeta);
            return (v_upper_orig->value(t, x) * b_der(one));
        };
        v_upper_ptr_ = std::make_shared<neumann_boundary_2d>(v_upper_trans);
    }

    // transform both horizontal:
    // horizontal lower:
    auto const h_lower_trans = [=](double t, double eta) -> double {
        auto const y = grid_2d::transformed_value_2(grid_transform_config, eta);
        return h_lower_orig->value(t, y);
    };

    if (auto hor_ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(h_upper_orig))
    {
        // horizontal upper:
        h_upper_trans = [=](double t, double eta) -> double {
            auto const y = grid_2d::transformed_value_2(grid_transform_config, eta);
            return (h_upper_orig->value(t, y) * a_der(one));
        };
        h_pair_ptr_ = std::make_pair(std::make_shared<dirichlet_boundary_2d>(h_lower_trans),
                                     std::make_shared<neumann_boundary_2d>(h_upper_trans));
    }
    else if (auto hor_ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(h_upper_orig))
    {
        // horizontal upper:
        h_upper_trans = [=](double t, double eta) -> double {
            auto const y = grid_2d::transformed_value_2(grid_transform_config, eta);
            return (h_upper_orig->value(t, y));
        };
        h_pair_ptr_ = std::make_pair(std::make_shared<dirichlet_boundary_2d>(h_lower_trans),
                                     std::make_shared<dirichlet_boundary_2d>(h_upper_trans));
    }
}

heston_boundary_transform::heston_boundary_transform(boundary_2d_ptr const &vertical_upper_boundary_ptr,
                                                     boundary_2d_pair const &horizontal_boundary_pair,
                                                     grid_transform_config_2d_ptr const grid_transform_config)
{
    initialize(vertical_upper_boundary_ptr, horizontal_boundary_pair, grid_transform_config);
}

heston_boundary_transform::~heston_boundary_transform()
{
}

boundary_2d_ptr const &heston_boundary_transform::vertical_upper() const
{
    return v_upper_ptr_;
}

boundary_2d_pair const &heston_boundary_transform::horizontal_pair() const
{
    return h_pair_ptr_;
}

} // namespace lss_pde_solvers
