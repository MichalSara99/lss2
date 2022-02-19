#if !defined(_ROBIN_3D_HPP_)
#define _ROBIN_3D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_robin_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using robin_3d_ptr = lss_boundary::robin_boundary_3d_ptr;
using robin_3d = lss_boundary::robin_boundary_3d;

struct robin_3d_builder
{
  private:
    std::function<double(double, double, double)> linear_value_;
    std::function<double(double, double, double)> value_;

  public:
    LSS_API explicit robin_3d_builder();

    LSS_API robin_3d_builder &linear_value(const std::function<double(double, double, double)> &linear_value);

    LSS_API robin_3d_builder &value(const std::function<double(double, double, double)> &value);

    LSS_API robin_3d_ptr build();
};

} // namespace lss

#endif ///_ROBIN_3D_HPP_
