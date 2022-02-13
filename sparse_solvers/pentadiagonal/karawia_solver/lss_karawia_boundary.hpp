/**

    @file      lss_karawia_boundary.hpp
    @brief     Boundary for Pentadiagonal Karawia Solver
    @details   ~
    @author    Michal Sara
    @date      13.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_KARAWIA_BOUNDARY_HPP_)
#define _LSS_KARAWIA_BOUNDARY_HPP_

#include <type_traits>
#include <vector>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../common/lss_enumerations.hpp"
#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"

namespace lss_karawia_solver
{

using lss_boundary::boundary_1d_pair;
using lss_boundary::boundary_2d_pair;
using lss_utility::sptr_t;

using sixtuple_t = std::tuple<double, double, double, double, double, double>;

class karawia_solver_boundary
{
  private:
    sixtuple_t lowest_sexta_;
    sixtuple_t lower_sexta_;
    sixtuple_t higher_sexta_;
    sixtuple_t highest_sexta_;
    double space_step_;
    double r0_, r1_, r_, rend_;
    std::size_t start_index_, end_index_;
    std::size_t discretization_size_;

    explicit karawia_solver_boundary() = delete;

    void initialise(boundary_1d_pair const &lowest_boundary, boundary_1d_pair const &lower_boundary, double time);

    void initialise(boundary_2d_pair const &lowest_boundary, boundary_2d_pair const &lower_boundary, double time,
                    double space_args);

    void finalise(boundary_1d_pair const &uppest_boundary, boundary_1d_pair const &upper_boundary, double time);

    void finalise(boundary_2d_pair const &uppest_boundary, boundary_2d_pair const &upper_boundary, double time,
                  double space_args);

  public:
    explicit karawia_solver_boundary(const std::size_t discretization_size, const double &space_step);

    ~karawia_solver_boundary();

    void set_lowest_sixtuple(const sixtuple_t &lowest_sexta);

    void set_lower_sixtuple(const sixtuple_t &lower_sexta);

    void set_highest_sixtuple(const sixtuple_t &highest_sexta);

    void set_higher_sixtuple(const sixtuple_t &higher_sexta);

    const std::tuple<double, double> init_coefficients(boundary_1d_pair const &lowest_boundary,
                                                       boundary_1d_pair const &lower_boundary, double time);

    const std::tuple<double, double> init_coefficients(boundary_2d_pair const &lowest_boundary,
                                                       boundary_2d_pair const &lower_boundary, double time,
                                                       double space_args);

    const std::tuple<double, double> final_coefficients(boundary_1d_pair const &uppest_boundary,
                                                        boundary_1d_pair const &upper_boundary, double time);

    const std::tuple<double, double> final_coefficients(boundary_2d_pair const &uppest_boundary,
                                                        boundary_2d_pair const &upper_boundary, double time,
                                                        double space_args);

    std::size_t start_index() const;

    std::size_t end_index() const;

    const double upper_boundary(boundary_1d_pair const &boundary, double time);

    const double upper_boundary(boundary_2d_pair const &boundary, double time, double space_args);

    const double lower_boundary(boundary_1d_pair const &boundary, double time);

    const double lower_boundary(boundary_2d_pair const &boundary, double time, double space_args);
};

using karawia_solver_boundary_ptr = sptr_t<karawia_solver_boundary>;

} // namespace lss_karawia_solver

#endif ///_LSS_KARAWIA_BOUNDARY_HPP_
