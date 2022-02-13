#if !defined(_NEUMANN_2D_HPP_)
#define _NEUMANN_2D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using neumann_2d_ptr = lss_boundary::neumann_boundary_2d_ptr;
using neumann_2d = lss_boundary::neumann_boundary_2d;

struct neumann_2d_builder
{
  private:
    std::function<double(double, double)> value_;

  public:
    LSS_API explicit neumann_2d_builder();

    LSS_API neumann_2d_builder &value(const std::function<double(double, double)> &value);

    LSS_API neumann_2d_ptr build();
};
} // namespace lss

#endif ///_NEUMANN_2D_HPP_
