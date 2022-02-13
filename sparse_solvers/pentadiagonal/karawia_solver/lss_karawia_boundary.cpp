#include "lss_karawia_boundary.hpp"

#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include <type_traits>
#include <vector>

namespace lss_karawia_solver
{

using lss_boundary::dirichlet_boundary_1d;
using lss_boundary::dirichlet_boundary_2d;

karawia_solver_boundary::karawia_solver_boundary(const std::size_t discretization_size, const double &space_step)
    : discretization_size_{discretization_size}, space_step_{space_step}
{
}

karawia_solver_boundary::~karawia_solver_boundary()
{
}

void karawia_solver_boundary::set_lowest_sixtuple(const sixtuple_t &lowest_sexta)
{
    lowest_sexta_ = lowest_sexta;
}

void karawia_solver_boundary::set_lower_sixtuple(const sixtuple_t &lower_sexta)
{
    lower_sexta_ = lower_sexta;
}

void karawia_solver_boundary::set_highest_sixtuple(const sixtuple_t &highest_sexta)
{
    highest_sexta_ = highest_sexta;
}

void karawia_solver_boundary::set_higher_sixtuple(const sixtuple_t &higher_sexta)
{
    higher_sexta_ = higher_sexta;
}

const std::tuple<double, double> karawia_solver_boundary::init_coefficients(boundary_1d_pair const &lowest_boundary,
                                                                            boundary_1d_pair const &lower_boundary,
                                                                            double time)
{
    initialise(lowest_boundary, lower_boundary, time);
    return std::make_tuple(r0_, r1_);
}

const std::tuple<double, double> karawia_solver_boundary::init_coefficients(boundary_2d_pair const &lowest_boundary,
                                                                            boundary_2d_pair const &lower_boundary,
                                                                            double time, double space_args)
{
    initialise(lowest_boundary, lower_boundary, time, space_args);
    return std::make_tuple(r0_, r1_);
}

const std::tuple<double, double> karawia_solver_boundary::final_coefficients(boundary_1d_pair const &uppest_boundary,
                                                                             boundary_1d_pair const &upper_boundary,
                                                                             double time)
{
    finalise(uppest_boundary, upper_boundary, time);
    return std::make_tuple(r_, rend_);
}

const std::tuple<double, double> karawia_solver_boundary::final_coefficients(boundary_2d_pair const &uppest_boundary,
                                                                             boundary_2d_pair const &upper_boundary,
                                                                             double time, double space_args)
{
    finalise(uppest_boundary, upper_boundary, time, space_args);
    return std::make_tuple(r_, rend_);
}

std::size_t karawia_solver_boundary::start_index() const
{
    return start_index_;
}

std::size_t karawia_solver_boundary::end_index() const
{
    return end_index_;
}

void karawia_solver_boundary::initialise(boundary_1d_pair const &lowest_boundary,
                                         boundary_1d_pair const &lower_boundary, double time)
{
    const auto a_2 = std::get<0>(lowest_sexta_);
    const auto b_2 = std::get<1>(lowest_sexta_);
    const auto f_2 = std::get<5>(lowest_sexta_);
    const auto a_3 = std::get<0>(lower_sexta_);
    const auto f_3 = std::get<5>(lower_sexta_);
    auto const lower_bnd = lower_boundary.first;
    auto const lowest_bnd = lowest_boundary.first;
    if (auto ptr_0 = std::dynamic_pointer_cast<dirichlet_boundary_1d>(lowest_bnd))
    {
        if (auto ptr_1 = std::dynamic_pointer_cast<dirichlet_boundary_1d>(lower_bnd))
        {
            const auto cst_val_0 = ptr_0->value(time);
            const auto cst_val_1 = ptr_1->value(time);
            start_index_ = 2;
            r0_ = f_2 - b_2 * cst_val_1 - a_2 * cst_val_0;
            r1_ = f_3 - a_3 * cst_val_1;
        }
        else
        {
            throw std::exception("Any other boundary type is not supported");
        }
    }
    else
    {
        throw std::exception("Any other boundary type is not supported");
    }
}

void karawia_solver_boundary::initialise(boundary_2d_pair const &lowest_boundary,
                                         boundary_2d_pair const &lower_boundary, double time, double space_args)
{
    const auto a_2 = std::get<0>(lowest_sexta_);
    const auto b_2 = std::get<1>(lowest_sexta_);
    const auto f_2 = std::get<5>(lowest_sexta_);
    const auto a_3 = std::get<0>(lower_sexta_);
    const auto f_3 = std::get<5>(lower_sexta_);
    auto const lower_bnd = lower_boundary.first;
    auto const lowest_bnd = lowest_boundary.first;
    if (auto ptr_0 = std::dynamic_pointer_cast<dirichlet_boundary_2d>(lowest_bnd))
    {
        if (auto ptr_1 = std::dynamic_pointer_cast<dirichlet_boundary_2d>(lower_bnd))
        {
            const auto cst_val_0 = ptr_0->value(time, space_args);
            const auto cst_val_1 = ptr_1->value(time, space_args);
            start_index_ = 2;
            r0_ = f_2 - b_2 * cst_val_1 - a_2 * cst_val_0;
            r1_ = f_3 - a_3 * cst_val_1;
        }
        else
        {
            throw std::exception("Any other boundary type is not supported");
        }
    }
    else
    {
        throw std::exception("Any other boundary type is not supported");
    }
}

void karawia_solver_boundary::finalise(boundary_1d_pair const &uppest_boundary, boundary_1d_pair const &upper_boundary,
                                       double time)
{
    const auto e = std::get<4>(higher_sexta_);
    const auto f = std::get<5>(higher_sexta_);
    const auto d_end = std::get<3>(highest_sexta_);
    const auto e_end = std::get<4>(highest_sexta_);
    const auto f_end = std::get<5>(highest_sexta_);
    auto const upper_bnd = upper_boundary.second;
    auto const uppest_bnd = uppest_boundary.second;
    if (auto ptr_end = std::dynamic_pointer_cast<dirichlet_boundary_1d>(uppest_bnd))
    {
        if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(upper_bnd))
        {
            const auto cst_val = ptr->value(time);
            const auto cst_val_end = ptr_end->value(time);
            end_index_ = discretization_size_ - 3;
            r_ = f - e * cst_val;
            rend_ = f_end - d_end * cst_val - e_end * cst_val_end;
        }
        else
        {
            throw std::exception("Any other boundary type is not supported");
        }
    }
    else
    {
        throw std::exception("Any other boundary type is not supported");
    }
}

void karawia_solver_boundary::finalise(boundary_2d_pair const &uppest_boundary, boundary_2d_pair const &upper_boundary,
                                       double time, double space_args)
{
    const auto e = std::get<4>(higher_sexta_);
    const auto f = std::get<5>(higher_sexta_);
    const auto d_end = std::get<3>(highest_sexta_);
    const auto e_end = std::get<4>(highest_sexta_);
    const auto f_end = std::get<5>(highest_sexta_);
    auto const upper_bnd = upper_boundary.second;
    auto const uppest_bnd = uppest_boundary.second;
    if (auto ptr_end = std::dynamic_pointer_cast<dirichlet_boundary_2d>(uppest_bnd))
    {
        if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(upper_bnd))
        {
            const auto cst_val = ptr->value(time, space_args);
            const auto cst_val_end = ptr_end->value(time, space_args);
            end_index_ = discretization_size_ - 3;
            r_ = f - e * cst_val;
            rend_ = f_end - d_end * cst_val - e_end * cst_val_end;
        }
        else
        {
            throw std::exception("Any other boundary type is not supported");
        }
    }
    else
    {
        throw std::exception("Any other boundary type is not supported");
    }
}

const double karawia_solver_boundary::upper_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(second_bnd))
    {
        ret = ptr->value(time);
    }

    return ret;
}

const double karawia_solver_boundary::upper_boundary(boundary_2d_pair const &boundary, double time, double space_args)
{
    double ret{};
    auto const second_bnd = boundary.second;
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(second_bnd))
    {
        ret = ptr->value(time, space_args);
    }

    return ret;
}

const double karawia_solver_boundary::lower_boundary(boundary_1d_pair const &boundary, double time)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_1d>(boundary.first))
    {
        ret = ptr->value(time);
    }
    return ret;
}

const double karawia_solver_boundary::lower_boundary(boundary_2d_pair const &boundary, double time, double space_args)
{
    double ret{};
    if (auto ptr = std::dynamic_pointer_cast<dirichlet_boundary_2d>(boundary.first))
    {
        ret = ptr->value(time, space_args);
    }
    return ret;
}

} // namespace lss_karawia_solver
