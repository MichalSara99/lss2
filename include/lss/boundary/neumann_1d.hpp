#if !defined(_NEUMANN_1D_HPP_)
#define _NEUMANN_1D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using neumann_1d_ptr = lss_boundary::neumann_boundary_1d_ptr;
using neumann_1d = lss_boundary::neumann_boundary_1d;

struct neumann_1d_builder
{
  private:
    double const_value_;
    std::function<double(double)> fun_value_;

  public:
    LSS_API explicit neumann_1d_builder();

    LSS_API neumann_1d_builder &value(double value);

    LSS_API neumann_1d_builder &value(const std::function<double(double)> &value);

    LSS_API neumann_1d_ptr build();
};
} // namespace lss

#endif ///_NEUMANN_1D_HPP_
