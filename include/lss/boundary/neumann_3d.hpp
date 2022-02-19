#if !defined(_NEUMANN_3D_HPP_)
#define _NEUMANN_3D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_neumann_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using neumann_3d_ptr = lss_boundary::neumann_boundary_3d_ptr;
using neumann_3d = lss_boundary::neumann_boundary_3d;

struct neumann_3d_builder
{
  private:
    std::function<double(double, double, double)> value_;

  public:
    LSS_API explicit neumann_3d_builder();

    LSS_API neumann_3d_builder &value(const std::function<double(double, double, double)> &value);

    LSS_API neumann_3d_ptr build();
};
} // namespace lss

#endif ///_NEUMANN_3D_HPP_
