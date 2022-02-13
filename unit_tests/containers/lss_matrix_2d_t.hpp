/**

    @file      lss_matrix_2d_t.hpp
    @brief     Unit tests for matrix_2d
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_MATRIX_2D_T_HPP_)
#define _MATRIX_2D_T_HPP_

#include "../../containers/lss_matrix_2d.hpp"
#include <iostream>
#include <random>
#include <string>

using lss_containers::matrix_2d;
namespace
{
void print_matrix_2d(matrix_2d const &m, std::string const &str = "")
{
    std::cout << str << ":\n";
    for (std::size_t r = 0; r < m.rows(); ++r)
    {
        for (std::size_t c = 0; c < m.columns(); ++c)
        {
            std::cout << m(r, c) << ",";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void print_array(std::valarray<double> const &a, std::string const &str = "")
{
    std::cout << str << ":\n";
    for (std::size_t c = 0; c < a.size(); ++c)
    {
        std::cout << a[c] << ",";
    }
    std::cout << "\n";
}

} // namespace

extern void basic_matrix_2d();

extern void slice_matrix_2d();

#endif ///_MATRIX_2D_T_HPP_
