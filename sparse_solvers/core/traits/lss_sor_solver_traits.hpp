/**

    @file      lss_sor_solver_traits.hpp
    @brief      Traits for SOR solvers
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_LSS_SOR_SOLVER_TRAITS_HPP_)
#define _LSS_SOR_SOLVER_TRAITS_HPP_

#include <typeinfo>

namespace lss_sor_solver_traits
{

template <typename fp_type> struct sor_solver_traits
{
};

template <> struct sor_solver_traits<double>
{
    static double tolerance()
    {
        return 1.e-18;
    }
    static std::size_t iteration_limit()
    {
        return 10'000;
    }
};

template <> struct sor_solver_traits<float>
{
    static float tolerance()
    {
        return 1.e-11f;
    }
    static std::size_t iteration_limit()
    {
        return 10'000;
    }
};

template <typename fp_type> struct sor_solver_cuda_traits
{
};

template <> struct sor_solver_cuda_traits<double>
{
    static double tolerance()
    {
        return 1.e-18;
    }
    static std::size_t iteration_limit()
    {
        return 100'000;
    }
};

template <> struct sor_solver_cuda_traits<float>
{
    static float tolerance()
    {
        return 1.e-11f;
    }
    static std::size_t iteration_limit()
    {
        return 100'000;
    }
};

} // namespace lss_sor_solver_traits

#endif ///_LSS_SOR_SOLVER_TRAITS_HPP_
