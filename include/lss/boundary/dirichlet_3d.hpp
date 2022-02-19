#if !defined(_DIRICHLET_2D_HPP_)
#define _DIRICHLET_2D_HPP_

#include <functional>

#include "../../../boundaries/lss_boundary.hpp"
#include "../../../boundaries/lss_dirichlet_boundary.hpp"
#include "../../../common/lss_macros.hpp"

namespace lss
{

using dirichlet_2d_ptr = lss_boundary::dirichlet_boundary_2d_ptr;
using dirichlet_2d = lss_boundary::dirichlet_boundary_2d;

struct dirichlet_2d_builder
{
  private:
    std::function<double(double, double)> value_;

  public:
    LSS_API explicit dirichlet_2d_builder();

    LSS_API dirichlet_2d_builder &value(const std::function<double(double, double)> &value);

    LSS_API dirichlet_2d_ptr build();
};
} // namespace lss

#endif ///_DIRICHLET_2D_HPP_
