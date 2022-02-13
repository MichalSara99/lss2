#include "lss_boundary_transform.hpp"

#include "../boundaries/lss_dirichlet_boundary.hpp"
#include "../boundaries/lss_neumann_boundary.hpp"
#include "../boundaries/lss_robin_boundary.hpp"
#include "../common/lss_macros.hpp"

namespace lss_transformation
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::robin_boundary_1d;

boundary_transform_1d::boundary_transform_1d(boundary_1d_pair const &boundary_pair,
                                             grid_transform_config_1d_ptr const grid_transform_config)
{
    initialize(boundary_pair, grid_transform_config);
}

boundary_transform_1d::~boundary_transform_1d()
{
}

boundary_1d_pair const &boundary_transform_1d::boundary_pair() const
{
    return pair_ptr_;
}

void boundary_transform_1d::initialize(boundary_1d_pair const &boundary_pair,
                                       grid_transform_config_1d_ptr const grid_transform_config)
{
    LSS_VERIFY(grid_transform_config, "grid_transform_config must not be null");
    auto const zero = static_cast<double>(0.0);
    auto const one = static_cast<double>(1.0);
    auto const &lower_orig = std::get<0>(boundary_pair);
    auto const &upper_orig = std::get<1>(boundary_pair);
    auto const &a_der = grid_transform_config->a_derivative();

    // first transform lower boundary:
    boundary_1d_ptr lower_ptr;
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(lower_orig))
    {
        // for Dirichlet no transform is necessary:
        lower_ptr = ptr;
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(lower_orig))
    {
        if (ptr->is_time_dependent())
        {
            auto const &lower_trans = [=](double t) -> double { return ptr->value(t) * a_der(zero); };
            lower_ptr = std::make_shared<neumann_boundary_1d>(lower_trans);
        }
        else
        {
            auto const lower_trans = ptr->value() * a_der(zero);
            lower_ptr = std::make_shared<neumann_boundary_1d>(lower_trans);
        }
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(lower_orig))
    {
        if (ptr->is_time_dependent())
        {
            auto const &lower_trans_lin = [=](double t) -> double { return ptr->linear_value(t) * a_der(zero); };
            auto const &lower_trans = [=](double t) -> double { return ptr->value(t) * a_der(zero); };
            lower_ptr = std::make_shared<robin_boundary_1d>(lower_trans_lin, lower_trans);
        }
        else
        {
            auto const lower_trans_lin = ptr->linear_value() * a_der(zero);
            auto const lower_trans = ptr->value() * a_der(zero);
            lower_ptr = std::make_shared<robin_boundary_1d>(lower_trans_lin, lower_trans);
        }
    }
    else
    {
        std::exception("Unreachable");
    }

    // second transform upper boundary:
    boundary_1d_ptr upper_ptr;
    if (auto const &ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(upper_orig))
    {
        // for Dirichlet no transform is necessary:
        upper_ptr = ptr;
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(upper_orig))
    {
        if (ptr->is_time_dependent())
        {
            auto const &upper_trans = [=](double t) -> double { return ptr->value(t) * a_der(one); };
            upper_ptr = std::make_shared<neumann_boundary_1d>(upper_trans);
        }
        else
        {
            auto const upper_trans = ptr->value() * a_der(one);
            upper_ptr = std::make_shared<neumann_boundary_1d>(upper_trans);
        }
    }
    else if (auto const &ptr = std::dynamic_pointer_cast<robin_boundary_1d>(upper_orig))
    {
        if (ptr->is_time_dependent())
        {
            auto const &upper_trans_lin = [=](double t) -> double { return ptr->linear_value(t) * a_der(one); };
            auto const &upper_trans = [=](double t) -> double { return ptr->value(t) * a_der(one); };
            upper_ptr = std::make_shared<robin_boundary_1d>(upper_trans_lin, upper_trans);
        }
        else
        {
            auto const upper_trans_lin = ptr->linear_value() * a_der(one);
            auto const upper_trans = ptr->value() * a_der(one);
            upper_ptr = std::make_shared<robin_boundary_1d>(upper_trans_lin, upper_trans);
        }
    }
    else
    {
        std::exception("Unreachable");
    }
    pair_ptr_ = std::make_pair(lower_ptr, upper_ptr);
}

} // namespace lss_transformation
