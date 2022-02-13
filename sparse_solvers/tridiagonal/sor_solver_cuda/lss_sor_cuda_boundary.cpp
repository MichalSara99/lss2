#include "lss_sor_cuda_boundary.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"

namespace lss_sor_solver_cuda
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

sor_cuda_boundary::sor_cuda_boundary(const std::size_t discretization_size, const double &space_step)
    : discretization_size_{discretization_size}, space_step_{space_step}
{
}

sor_cuda_boundary::~sor_cuda_boundary()
{
}

void sor_cuda_boundary::set_lowest_quad(const quad_t &lowest_quad)
{
    lowest_quad_ = lowest_quad;
}

void sor_cuda_boundary::set_lower_quad(const quad_t &lower_quad)
{
    lower_quad_ = lower_quad;
}

void sor_cuda_boundary::set_highest_quad(const quad_t &highest_quad)
{
    highest_quad_ = highest_quad;
}

void sor_cuda_boundary::set_higher_quad(const quad_t &higher_quad)
{
    higher_quad_ = higher_quad;
}

const triplet_t sor_cuda_boundary::init_coefficients(boundary_1d_pair const &boundary, double time)
{
    initialise(boundary, time);
    return std::make_tuple(b_init_, c_init_, f_init_);
}

const triplet_t sor_cuda_boundary::init_coefficients(boundary_2d_pair const &boundary, double time, double space_arg)
{
    initialise(boundary, time, space_arg);
    return std::make_tuple(b_init_, c_init_, f_init_);
}

const triplet_t sor_cuda_boundary::init_coefficients(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                                     double space_2_arg)
{
    initialise(boundary, time, space_1_arg, space_2_arg);
    return std::make_tuple(b_init_, c_init_, f_init_);
}

const triplet_t sor_cuda_boundary::final_coefficients(boundary_1d_pair const &boundary, double time)
{
    finalise(boundary, time);
    return std::make_tuple(a_end_, b_end_, f_end_);
}

const triplet_t sor_cuda_boundary::final_coefficients(boundary_2d_pair const &boundary, double time, double space_arg)
{
    finalise(boundary, time, space_arg);
    return std::make_tuple(a_end_, b_end_, f_end_);
}

const triplet_t sor_cuda_boundary::final_coefficients(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                                      double space_2_arg)
{
    finalise(boundary, time, space_1_arg, space_2_arg);
    return std::make_tuple(a_end_, b_end_, f_end_);
}

std::size_t sor_cuda_boundary::start_index() const
{
    return start_index_;
}

std::size_t sor_cuda_boundary::end_index() const
{
    return end_index_;
}

void sor_cuda_boundary::initialise(boundary_1d_pair const &boundary, double time)
{
    const auto a_0 = std::get<0>(lowest_quad_);
    const auto b_0 = std::get<1>(lowest_quad_);
    const auto c_0 = std::get<2>(lowest_quad_);
    const auto f_0 = std::get<3>(lowest_quad_);
    const auto a_1 = std::get<0>(lower_quad_);
    const auto b_1 = std::get<1>(lower_quad_);
    const auto c_1 = std::get<2>(lower_quad_);
    const auto f_1 = std::get<3>(lower_quad_);
    const double two = static_cast<double>(2.0);

    auto const &first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(first_bnd))
    {
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 1;
        b_init_ = b_1;
        c_init_ = c_1;
        f_init_ = f_1 - a_1 * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 0;
        b_init_ = b_0;
        c_init_ = a_0 + c_0;
        f_init_ = f_0 - a_0 * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const auto lin_val =
            two * space_step_ * (ptr->is_time_dependent() ? ptr->linear_value(time) : ptr->linear_value());
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 0;
        b_init_ = b_0 + a_0 * cst_val;
        c_init_ = a_0 + c_0;
        f_init_ = f_0 - a_0 * cst_val;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void sor_cuda_boundary::initialise(boundary_2d_pair const &boundary, double time, double space_arg)
{
    const auto a_0 = std::get<0>(lowest_quad_);
    const auto b_0 = std::get<1>(lowest_quad_);
    const auto c_0 = std::get<2>(lowest_quad_);
    const auto f_0 = std::get<3>(lowest_quad_);
    const auto a_1 = std::get<0>(lower_quad_);
    const auto b_1 = std::get<1>(lower_quad_);
    const auto c_1 = std::get<2>(lower_quad_);
    const auto f_1 = std::get<3>(lower_quad_);
    const double two = static_cast<double>(2.0);

    auto const &first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_arg);
        start_index_ = 1;
        b_init_ = b_1;
        c_init_ = c_1;
        f_init_ = f_1 - a_1 * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(first_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        start_index_ = 0;
        b_init_ = b_0;
        c_init_ = a_0 + c_0;
        f_init_ = f_0 - a_0 * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_2d>(first_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        start_index_ = 0;
        b_init_ = b_0 + a_0 * cst_val;
        c_init_ = a_0 + c_0;
        f_init_ = f_0 - a_0 * cst_val;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void sor_cuda_boundary::initialise(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                   double space_2_arg)
{
    const auto a_0 = std::get<0>(lowest_quad_);
    const auto b_0 = std::get<1>(lowest_quad_);
    const auto c_0 = std::get<2>(lowest_quad_);
    const auto f_0 = std::get<3>(lowest_quad_);
    const auto a_1 = std::get<0>(lower_quad_);
    const auto b_1 = std::get<1>(lower_quad_);
    const auto c_1 = std::get<2>(lower_quad_);
    const auto f_1 = std::get<3>(lower_quad_);
    const double two = static_cast<double>(2.0);

    auto const &first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 1;
        b_init_ = b_1;
        c_init_ = c_1;
        f_init_ = f_1 - a_1 * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(first_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 0;
        b_init_ = b_0;
        c_init_ = a_0 + c_0;
        f_init_ = f_0 - a_0 * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_3d>(first_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_1_arg, space_2_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 0;
        b_init_ = b_0 + a_0 * cst_val;
        c_init_ = a_0 + c_0;
        f_init_ = f_0 - a_0 * cst_val;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void sor_cuda_boundary::finalise(boundary_1d_pair const &boundary, double time)
{
    const auto a = std::get<0>(higher_quad_);
    const auto b = std::get<1>(higher_quad_);
    const auto c = std::get<2>(higher_quad_);
    const auto f = std::get<3>(higher_quad_);
    const auto a_end = std::get<0>(highest_quad_);
    const auto b_end = std::get<1>(highest_quad_);
    const auto c_end = std::get<2>(highest_quad_);
    const auto f_end = std::get<3>(highest_quad_);
    const double two = static_cast<double>(2.0);

    auto const &second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        end_index_ = discretization_size_ - 2;
        a_end_ = a;
        b_end_ = b;
        f_end_ = f - c * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        end_index_ = discretization_size_ - 1;
        a_end_ = a_end + c_end;
        b_end_ = b_end;
        f_end_ = f_end + c_end * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const auto lin_val =
            two * space_step_ * (ptr->is_time_dependent() ? ptr->linear_value(time) : ptr->linear_value());
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        end_index_ = discretization_size_ - 1;
        a_end_ = a_end + c_end;
        b_end_ = b_end - c_end * lin_val;
        f_end_ = f_end + c_end * cst_val;
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

void sor_cuda_boundary::finalise(boundary_2d_pair const &boundary, double time, double space_arg)
{
    const auto a = std::get<0>(higher_quad_);
    const auto b = std::get<1>(higher_quad_);
    const auto c = std::get<2>(higher_quad_);
    const auto f = std::get<3>(higher_quad_);
    const auto a_end = std::get<0>(highest_quad_);
    const auto b_end = std::get<1>(highest_quad_);
    const auto c_end = std::get<2>(highest_quad_);
    const auto f_end = std::get<3>(highest_quad_);
    const double two = static_cast<double>(2.0);

    auto const &second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        const auto cst_val = ptr->value(time, space_arg);
        end_index_ = discretization_size_ - 2;
        a_end_ = a;
        b_end_ = b;
        f_end_ = f - c * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        end_index_ = discretization_size_ - 1;
        a_end_ = a_end + c_end;
        b_end_ = b_end;
        f_end_ = f_end + c_end * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_2d>(second_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        end_index_ = discretization_size_ - 1;
        a_end_ = a_end + c_end;
        b_end_ = b_end - c_end * lin_val;
        f_end_ = f_end + c_end * cst_val;
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

void sor_cuda_boundary::finalise(boundary_3d_pair const &boundary, double time, double space_1_arg, double space_2_arg)
{
    const auto a = std::get<0>(higher_quad_);
    const auto b = std::get<1>(higher_quad_);
    const auto c = std::get<2>(higher_quad_);
    const auto f = std::get<3>(higher_quad_);
    const auto a_end = std::get<0>(highest_quad_);
    const auto b_end = std::get<1>(highest_quad_);
    const auto c_end = std::get<2>(highest_quad_);
    const auto f_end = std::get<3>(highest_quad_);
    const double two = static_cast<double>(2.0);

    auto const &second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(second_bnd))
    {
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        end_index_ = discretization_size_ - 2;
        a_end_ = a;
        b_end_ = b;
        f_end_ = f - c * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        end_index_ = discretization_size_ - 1;
        a_end_ = a_end + c_end;
        b_end_ = b_end;
        f_end_ = f_end + c_end * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_3d>(second_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_1_arg, space_2_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        end_index_ = discretization_size_ - 1;
        a_end_ = a_end + c_end;
        b_end_ = b_end - c_end * lin_val;
        f_end_ = f_end + c_end * cst_val;
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

const double sor_cuda_boundary::upper_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(boundary.second))
    {
        ret = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
    }
    return ret;
}

const double sor_cuda_boundary::upper_boundary(boundary_2d_pair const &boundary, double time, double space_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(boundary.second))
    {
        ret = ptr->value(time, space_arg);
    }
    return ret;
}

const double sor_cuda_boundary::upper_boundary(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                               double space_2_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(boundary.second))
    {
        ret = ptr->value(time, space_1_arg, space_2_arg);
    }
    return ret;
}

const double sor_cuda_boundary::lower_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(boundary.first))
    {
        ret = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
    }
    return ret;
}

const double sor_cuda_boundary::lower_boundary(boundary_2d_pair const &boundary, double time, double space_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(boundary.first))
    {
        ret = ptr->value(time, space_arg);
    }
    return ret;
}

const double sor_cuda_boundary::lower_boundary(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                               double space_2_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(boundary.first))
    {
        ret = ptr->value(time, space_1_arg, space_2_arg);
    }
    return ret;
}

} // namespace lss_sor_solver_cuda
