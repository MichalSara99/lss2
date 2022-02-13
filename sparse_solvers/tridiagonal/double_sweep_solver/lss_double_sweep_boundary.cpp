#include "lss_double_sweep_boundary.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss_double_sweep_solver
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::dirichlet_boundary_2d;
using lss_boundary::dirichlet_boundary_3d;
using lss_boundary::neumann_boundary_1d;
using lss_boundary::neumann_boundary_2d;
using lss_boundary::neumann_boundary_3d;
using lss_boundary::robin_boundary_1d;
using lss_boundary::robin_boundary_2d;
using lss_boundary::robin_boundary_3d;

double_sweep_boundary::double_sweep_boundary(const std::size_t &discretization_size, const double &space_step)
    : discretization_size_{discretization_size}, space_step_{space_step}
{
}

double_sweep_boundary::~double_sweep_boundary()
{
}

void double_sweep_boundary::set_low_quad(const quad_t &quad)
{
    low_quad_ = quad;
}

std::size_t double_sweep_boundary::start_index() const
{
    return start_index_;
}

std::size_t double_sweep_boundary::end_index(boundary_1d_pair const &boundary) const
{
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(boundary.second))
    {
        return (discretization_size_ - 2);
    }
    return (discretization_size_ - 1);
}

std::size_t double_sweep_boundary::end_index(boundary_2d_pair const &boundary) const
{
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(boundary.second))
    {
        return (discretization_size_ - 2);
    }
    return (discretization_size_ - 1);
}

std::size_t double_sweep_boundary::end_index(boundary_3d_pair const &boundary) const
{
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(boundary.second))
    {
        return (discretization_size_ - 2);
    }
    return (discretization_size_ - 1);
}

void double_sweep_boundary::initialise(boundary_1d_pair const &boundary, double time)
{
    const auto a = std::get<0>(low_quad_);
    const auto b = std::get<1>(low_quad_);
    const auto c = std::get<2>(low_quad_);
    const auto f = std::get<3>(low_quad_);
    const double two = static_cast<double>(2.0);
    const double mone = static_cast<double>(-1.0);
    auto const &first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(first_bnd))
    {
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 1;
        k_ = cst_val;
        l_ = double{};
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 0;
        k_ = (f - a * space_step_ * two * cst_val) / b;
        l_ = mone * (a + c) / b;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const auto lin_val = (ptr->is_time_dependent() ? ptr->linear_value(time) : ptr->linear_value());
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 0;
        const auto tmp = b + a * space_step_ * two * lin_val;
        k_ = (f - a * space_step_ * two * cst_val) / tmp;
        l_ = mone * (a + c) / tmp;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void double_sweep_boundary::initialise(boundary_2d_pair const &boundary, double time, double space_arg)
{
    const auto a = std::get<0>(low_quad_);
    const auto b = std::get<1>(low_quad_);
    const auto c = std::get<2>(low_quad_);
    const auto f = std::get<3>(low_quad_);
    const double two = static_cast<double>(2.0);
    const double mone = static_cast<double>(-1.0);
    auto const &first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_arg);
        start_index_ = 1;
        k_ = cst_val;
        l_ = float{};
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_arg);
        start_index_ = 0;
        k_ = (f - a * space_step_ * two * cst_val) / b;
        l_ = mone * (a + c) / b;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_2d>(first_bnd))
    {
        const auto lin_val = ptr->linear_value(time, space_arg);
        const auto cst_val = ptr->value(time, space_arg);
        start_index_ = 0;
        const auto tmp = b + a * space_step_ * two * lin_val;
        k_ = (f - a * space_step_ * two * cst_val) / tmp;
        l_ = mone * (a + c) / tmp;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void double_sweep_boundary::initialise(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                       double space_2_arg)
{
    const auto a = std::get<0>(low_quad_);
    const auto b = std::get<1>(low_quad_);
    const auto c = std::get<2>(low_quad_);
    const auto f = std::get<3>(low_quad_);
    const double two = static_cast<double>(2.0);
    const double mone = static_cast<double>(-1.0);
    auto const &first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 1;
        k_ = cst_val;
        l_ = float{};
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 0;
        k_ = (f - a * space_step_ * two * cst_val) / b;
        l_ = mone * (a + c) / b;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_3d>(first_bnd))
    {
        const auto lin_val = ptr->linear_value(time, space_1_arg, space_2_arg);
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 0;
        const auto tmp = b + a * space_step_ * two * lin_val;
        k_ = (f - a * space_step_ * two * cst_val) / tmp;
        l_ = mone * (a + c) / tmp;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void double_sweep_boundary::finalise(boundary_1d_pair const &boundary, const double &k_nm1, const double &k_n,
                                     const double &l_nm1, const double &l_n, double time)
{
    const double two = static_cast<double>(2.0);
    const double one = static_cast<double>(1.0);
    auto const &second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        upper_ = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        upper_ = (l_n * (k_nm1 - cst_val) + k_n) / (one - l_n * l_nm1);
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const auto lin_val =
            two * space_step_ * (ptr->is_time_dependent() ? ptr->linear_value(time) : ptr->linear_value());
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        upper_ = (l_n * (k_nm1 - cst_val) + k_n) / (one - l_n * (l_nm1 - lin_val));
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

void double_sweep_boundary::finalise(boundary_2d_pair const &boundary, const double &k_nm1, const double &k_n,
                                     const double &l_nm1, const double &l_n, double time, double space_args)
{
    const double two = static_cast<double>(2.0);
    const double one = static_cast<double>(1.0);
    auto const &second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        upper_ = ptr->value(time, space_args);
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_args);
        upper_ = (l_n * (k_nm1 - cst_val) + k_n) / (one - l_n * l_nm1);
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_2d>(second_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_args);
        const auto cst_val = two * space_step_ * ptr->value(time, space_args);
        upper_ = (l_n * (k_nm1 - cst_val) + k_n) / (one - l_n * (l_nm1 - lin_val));
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

void double_sweep_boundary::finalise(boundary_3d_pair const &boundary, const double &k_nm1, const double &k_n,
                                     const double &l_nm1, const double &l_n, double time, double space_1_arg,
                                     double space_2_arg)
{
    const double two = static_cast<double>(2.0);
    const double one = static_cast<double>(1.0);
    auto const &second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(second_bnd))
    {
        upper_ = ptr->value(time, space_1_arg, space_2_arg);
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        upper_ = (l_n * (k_nm1 - cst_val) + k_n) / (one - l_n * l_nm1);
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_3d>(second_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_1_arg, space_2_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        upper_ = (l_n * (k_nm1 - cst_val) + k_n) / (one - l_n * (l_nm1 - lin_val));
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

const std::pair<double, double> double_sweep_boundary::coefficients(boundary_1d_pair const &boundary, double time)
{
    initialise(boundary, time);
    return std::make_pair(k_, l_);
}

const std::pair<double, double> double_sweep_boundary::coefficients(boundary_2d_pair const &boundary, double time,
                                                                    double space_arg)
{
    initialise(boundary, time, space_arg);
    return std::make_pair(k_, l_);
}

const std::pair<double, double> double_sweep_boundary::coefficients(boundary_3d_pair const &boundary, double time,
                                                                    double space_1_arg, double space_2_arg)
{
    initialise(boundary, time, space_1_arg, space_2_arg);
    return std::make_pair(k_, l_);
}

const double double_sweep_boundary::upper_boundary(boundary_1d_pair const &boundary, const double &k_nm1,
                                                   const double &k_n, const double &l_nm1, const double &l_n,
                                                   double time)
{
    finalise(boundary, k_nm1, k_n, l_nm1, l_n, time);
    return upper_;
}

const double double_sweep_boundary::upper_boundary(boundary_2d_pair const &boundary, const double &k_nm1,
                                                   const double &k_n, const double &l_nm1, const double &l_n,
                                                   double time, double space_arg)
{
    finalise(boundary, k_nm1, k_n, l_nm1, l_n, time, space_arg);
    return upper_;
}

const double double_sweep_boundary::upper_boundary(boundary_3d_pair const &boundary, const double &k_nm1,
                                                   const double &k_n, const double &l_nm1, const double &l_n,
                                                   double time, double space_1_arg, double space_2_arg)
{
    finalise(boundary, k_nm1, k_n, l_nm1, l_n, time, space_1_arg, space_2_arg);
    return upper_;
}

const double double_sweep_boundary::lower_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(boundary.first))
    {
        ret = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
    }
    return ret;
}

const double double_sweep_boundary::lower_boundary(boundary_2d_pair const &boundary, double time, double space_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(boundary.first))
    {
        ret = ptr->value(time, space_arg);
    }
    return ret;
}

const double double_sweep_boundary::lower_boundary(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                                   double space_2_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(boundary.first))
    {
        ret = ptr->value(time, space_1_arg, space_2_arg);
    }
    return ret;
}

// double precision specialization:

} // namespace lss_double_sweep_solver
