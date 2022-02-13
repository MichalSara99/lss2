/**

    @file      lss_enumerations.hpp
    @brief     All library enumerations
    @details   ~
    @author    Michal Sara
    @date      14.12.2021
    @copyright © Michal Sara, 2021. All right reserved.

**/
#pragma once
#if !defined(_LSS_ENUMERATIONS_HPP_)
#define _LSS_ENUMERATIONS_HPP_

#include "lss_macros.hpp"

namespace lss_enumerations
{

/**
    @enum  lss_enumerations::factorization_enum
    @brief Factorization for CUDA solvers
**/
enum class factorization_enum
{
    QRMethod,
    LUMethod,
    CholeskyMethod,
    None,
};

/**
    @enum  lss_enumerations::splitting_method_enum
    @brief Methods for 2D PDE solvers
**/
enum class splitting_method_enum
{
    DouglasRachford,
    CraigSneyd,
    ModifiedCraigSneyd,
    HundsdorferVerwer,
};

/**
    @enum  lss_enumerations::tridiagonal_method_enum
    @brief Tridiagonal solvers
**/
enum class tridiagonal_method_enum
{
    CUDASolver,
    DoubleSweepSolver,
    SORSolver,
    ThomasLUSolver,
};

/**
    @enum  lss_enumerations::dimension_enum
    @brief Dimension enum
**/
enum class dimension_enum
{
    One,
    Two,
};

/**
    @enum  lss_enumerations::grid_enum
    @brief Spacial spacing for solvers
**/
enum class grid_enum
{
    Uniform,
    Nonuniform,
};

/**
    @enum  lss_enumerations::traverse_direction_enum
    @brief Traverse disrection for PDE solvers
**/
enum class traverse_direction_enum
{
    Forward,
    Backward,
};

/**
    @enum  lss_enumerations::memory_space_enum
    @brief Memory space for solvers
**/
enum class memory_space_enum
{
    Host,
    Device
};

/**
    @enum  lss_enumerations::flat_matrix_sort_enum
    @brief Flat matrix sorting
**/
enum class flat_matrix_sort_enum
{
    RowMajor,
    ColumnMajor
};

/**
    @enum  lss_enumerations::implicit_pde_schemes_enum
    @brief Implicit 1D, 2D and 3D PDE schemes
**/
enum class implicit_pde_schemes_enum
{
    Euler,
    CrankNicolson,
};

/**
    @enum  lss_enumerations::explicit_pde_schemes_enum
    @brief Explicit 1D PDE schemes
**/
enum class explicit_pde_schemes_enum
{
    Euler,
    ADEBarakatClark,
    ADESaulyev
};

} // namespace lss_enumerations

#endif ///_LSS_ENUMERATIONS_HPP_
