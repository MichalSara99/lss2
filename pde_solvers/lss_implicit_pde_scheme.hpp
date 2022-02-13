/**

    @file      lss_implicit_pde_scheme.hpp
    @brief     Implicit PDE scheme structure
    @details   ~
    @author    Michal Sara
    @date      30.01.2022
    @copyright © Michal Sara, 2022. All right reserved.

**/
#pragma once

#if !defined(_LSS_IMPLICIT_PDE_SCHEME_HPP_)
#define _LSS_IMPLICIT_PDE_SCHEME_HPP_

#include "../common/lss_enumerations.hpp"
#include "../common/lss_utility.hpp"

namespace lss_pde_solvers
{

using lss_enumerations::implicit_pde_schemes_enum;
using lss_utility::sptr_t;

/**
    @struct implicit_pde_scheme
    @brief structure representing implicit PDE scheme
**/
struct implicit_pde_scheme
{
  private:
    double value_;

    explicit implicit_pde_scheme();

  public:
    explicit implicit_pde_scheme(double value);

    explicit implicit_pde_scheme(implicit_pde_schemes_enum implicit_pde_enum);

    virtual ~implicit_pde_scheme();

    /**
        @brief Value for the implicit scheme
        @retval theta value for the pde scheme
    **/
    LSS_API double value() const;
};

using implicit_pde_scheme_ptr = sptr_t<implicit_pde_scheme>;

} // namespace lss_pde_solvers

#endif ///_LSS_IMPLICIT_PDE_SCHEME_HPP_
