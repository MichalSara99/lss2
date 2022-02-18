/**

    @file      lss_matrix_3d_t.hpp
    @brief     Unit tests for matrix_3d
    @details   ~
    @author    Michal Sara
    @date      12.02.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once
#if !defined(_MATRIX_3D_T_HPP_)
#define _MATRIX_3D_T_HPP_

#include "../../containers/lss_matrix_3d.hpp"
#include <iostream>
#include <random>
#include <string>

using lss_containers::matrix_3d;

namespace
{
void print_matrix_3d(matrix_3d const &m, std::string const &str = "")
{
    std::cout << str << ":\n";
    for (std::size_t l = 0; l < m.layers(); ++l)
    {
        for (std::size_t r = 0; r < m.rows(); ++r)
        {
            for (std::size_t c = 0; c < m.columns(); ++c)
            {
                std::cout << m(r, c, l) << ",";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

} // namespace

extern void basic_matrix_3d();

extern void slice_matrix_3d();

extern void slice_planes_matrix_3d();

#endif ///_MATRIX_3D_T_HPP_
