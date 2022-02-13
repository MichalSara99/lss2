#include "lss_thomas_lu_boundary.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"

namespace lss_thomas_lu_solver

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

thomas_lu_solver_boundary::thomas_lu_solver_boundary(const std::size_t discretization_size, const double &space_step)
    : discretization_size_{discretization_size}, space_step_{space_step}
{
}

thomas_lu_solver_boundary::~thomas_lu_solver_boundary()
{
}

void thomas_lu_solver_boundary::set_lowest_quad(const quad_t &lowest_quad)
{
    lowest_quad_ = lowest_quad;
}

void thomas_lu_solver_boundary::set_lower_quad(const quad_t &lower_quad)
{
    lower_quad_ = lower_quad;
}

void thomas_lu_solver_boundary::set_highest_quad(const quad_t &highest_quad)
{
    highest_quad_ = highest_quad;
}

void thomas_lu_solver_boundary::set_higher_quad(const quad_t &higher_quad)
{
    higher_quad_ = higher_quad;
}

const quad_t thomas_lu_solver_boundary::init_coefficients(boundary_1d_pair const &boundary, double time)
{
    initialise(boundary, time);
    return std::make_tuple(beta_, gamma_, r_, z_);
}

const quad_t thomas_lu_solver_boundary::init_coefficients(boundary_2d_pair const &boundary, double time,
                                                          double space_arg)
{
    initialise(boundary, time, space_arg);
    return std::make_tuple(beta_, gamma_, r_, z_);
}

const quad_t thomas_lu_solver_boundary::init_coefficients(boundary_3d_pair const &boundary, double time,
                                                          double space_1_arg, double space_2_arg)
{
    initialise(boundary, time, space_1_arg, space_2_arg);
    return std::make_tuple(beta_, gamma_, r_, z_);
}

const triplet_t thomas_lu_solver_boundary::final_coefficients(boundary_1d_pair const &boundary, double time)
{
    finalise(boundary, time);
    return std::make_tuple(alpha_n_, beta_n_, r_n_);
}

const triplet_t thomas_lu_solver_boundary::final_coefficients(boundary_2d_pair const &boundary, double time,
                                                              double space_arg)
{
    finalise(boundary, time, space_arg);
    return std::make_tuple(alpha_n_, beta_n_, r_n_);
}

const triplet_t thomas_lu_solver_boundary::final_coefficients(boundary_3d_pair const &boundary, double time,
                                                              double space_1_arg, double space_2_arg)
{
    finalise(boundary, time, space_1_arg, space_2_arg);
    return std::make_tuple(alpha_n_, beta_n_, r_n_);
}

std::size_t thomas_lu_solver_boundary::start_index() const
{
    return start_index_;
}

std::size_t thomas_lu_solver_boundary::end_index() const
{
    return end_index_;
}

void thomas_lu_solver_boundary::initialise(boundary_1d_pair const &boundary, double time)
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
    auto const first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(first_bnd))
    {
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 1;
        beta_ = b_1;
        gamma_ = c_1 / beta_;
        r_ = f_1 - a_1 * cst_val;
        z_ = r_ / beta_;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(first_bnd))
    {
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 0;
        beta_ = b_0;
        gamma_ = (a_0 + c_0) / beta_;
        r_ = f_0 - a_0 * cst_val;
        z_ = r_ / beta_;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_1d>(first_bnd))
    {
        const auto lin_val =
            two * space_step_ * (ptr->is_time_dependent() ? ptr->linear_value(time) : ptr->linear_value());
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        start_index_ = 0;
        beta_ = b_0 + a_0 * lin_val;
        gamma_ = (a_0 + c_0) / beta_;
        r_ = f_0 - a_0 * cst_val;
        z_ = r_ / beta_;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void thomas_lu_solver_boundary::initialise(boundary_2d_pair const &boundary, double time, double space_arg)
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
    auto const first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_arg);
        start_index_ = 1;
        beta_ = b_1;
        gamma_ = c_1 / beta_;
        r_ = f_1 - a_1 * cst_val;
        z_ = r_ / beta_;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(first_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        start_index_ = 0;
        beta_ = b_0;
        gamma_ = (a_0 + c_0) / beta_;
        r_ = f_0 - a_0 * cst_val;
        z_ = r_ / beta_;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_2d>(first_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        start_index_ = 0;
        beta_ = b_0 + a_0 * lin_val;
        gamma_ = (a_0 + c_0) / beta_;
        r_ = f_0 - a_0 * cst_val;
        z_ = r_ / beta_;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void thomas_lu_solver_boundary::initialise(boundary_3d_pair const &boundary, double time, double space_1_arg,
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
    auto const first_bnd = boundary.first;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(first_bnd))
    {
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 1;
        beta_ = b_1;
        gamma_ = c_1 / beta_;
        r_ = f_1 - a_1 * cst_val;
        z_ = r_ / beta_;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(first_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 0;
        beta_ = b_0;
        gamma_ = (a_0 + c_0) / beta_;
        r_ = f_0 - a_0 * cst_val;
        z_ = r_ / beta_;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_3d>(first_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_1_arg, space_2_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        start_index_ = 0;
        beta_ = b_0 + a_0 * lin_val;
        gamma_ = (a_0 + c_0) / beta_;
        r_ = f_0 - a_0 * cst_val;
        z_ = r_ / beta_;
    }
    else
    {
        // throw here unrecognized lower boundary
    }
}

void thomas_lu_solver_boundary::finalise(boundary_1d_pair const &boundary, double time)
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
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        const auto cst_val = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        end_index_ = discretization_size_ - 2;
        alpha_n_ = a;
        beta_n_ = b;
        r_n_ = f - c * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_1d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        end_index_ = discretization_size_ - 1;
        alpha_n_ = a_end + c_end;
        beta_n_ = b_end;
        r_n_ = f_end + c_end * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_1d>(second_bnd))
    {
        const auto lin_val =
            two * space_step_ * (ptr->is_time_dependent() ? ptr->linear_value(time) : ptr->linear_value());
        const auto cst_val = two * space_step_ * (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
        end_index_ = discretization_size_ - 1;
        alpha_n_ = a_end + c_end;
        beta_n_ = b_end - c_end * lin_val;
        r_n_ = f_end + c_end * cst_val;
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

void thomas_lu_solver_boundary::finalise(boundary_2d_pair const &boundary, double time, double space_arg)
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
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        const auto cst_val = ptr->value(time, space_arg);
        end_index_ = discretization_size_ - 2;
        alpha_n_ = a;
        beta_n_ = b;
        r_n_ = f - c * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_2d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        end_index_ = discretization_size_ - 1;
        alpha_n_ = a_end + c_end;
        beta_n_ = b_end;
        r_n_ = f_end + c_end * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_2d>(second_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_arg);
        end_index_ = discretization_size_ - 1;
        alpha_n_ = a_end + c_end;
        beta_n_ = b_end - c_end * lin_val;
        r_n_ = f_end + c_end * cst_val;
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

void thomas_lu_solver_boundary::finalise(boundary_3d_pair const &boundary, double time, double space_1_arg,
                                         double space_2_arg)
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
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(second_bnd))
    {
        const auto cst_val = ptr->value(time, space_1_arg, space_2_arg);
        end_index_ = discretization_size_ - 2;
        alpha_n_ = a;
        beta_n_ = b;
        r_n_ = f - c * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<neumann_boundary_3d>(second_bnd))
    {
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        end_index_ = discretization_size_ - 1;
        alpha_n_ = a_end + c_end;
        beta_n_ = b_end;
        r_n_ = f_end + c_end * cst_val;
    }
    else if (auto ptr = std::dynamic_pointer_cast<robin_boundary_3d>(second_bnd))
    {
        const auto lin_val = two * space_step_ * ptr->linear_value(time, space_1_arg, space_2_arg);
        const auto cst_val = two * space_step_ * ptr->value(time, space_1_arg, space_2_arg);
        end_index_ = discretization_size_ - 1;
        alpha_n_ = a_end + c_end;
        beta_n_ = b_end - c_end * lin_val;
        r_n_ = f_end + c_end * cst_val;
    }
    else
    {
        // throw here unrecognized upper boundary
    }
}

const double thomas_lu_solver_boundary::upper_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        ret = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
    }

    return ret;
}

const double thomas_lu_solver_boundary::upper_boundary(boundary_2d_pair const &boundary, double time, double space_arg)
{
    double ret{};
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        ret = ptr->value(time, space_arg);
    }

    return ret;
}

const double thomas_lu_solver_boundary::upper_boundary(boundary_3d_pair const &boundary, double time,
                                                       double space_1_arg, double space_2_arg)
{
    double ret{};
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(second_bnd))
    {
        ret = ptr->value(time, space_1_arg, space_2_arg);
    }

    return ret;
}

const double thomas_lu_solver_boundary::lower_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(boundary.first))
    {
        ret = (ptr->is_time_dependent() ? ptr->value(time) : ptr->value());
    }
    return ret;
}

const double thomas_lu_solver_boundary::lower_boundary(boundary_2d_pair const &boundary, double time, double space_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(boundary.first))
    {
        ret = ptr->value(time, space_arg);
    }
    return ret;
}

const double thomas_lu_solver_boundary::lower_boundary(boundary_3d_pair const &boundary, double time,
                                                       double space_1_arg, double space_2_arg)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_3d>(boundary.first))
    {
        ret = ptr->value(time, space_1_arg, space_2_arg);
    }
    return ret;
}

} // namespace lss_thomas_lu_solver
