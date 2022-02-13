#if !defined(_CONTAINER_3D_HPP_)
#define _CONTAINER_3D_HPP_

#include <functional>

#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_container_3d.hpp"

namespace lss
{

using rmatrix_3d = lss_containers::rmatrix_3d;
using rmatrix_3d_ptr = lss_containers::rmatrix_3d_ptr;

using cmatrix_3d = lss_containers::cmatrix_3d;
using cmatrix_3d_ptr = lss_containers::cmatrix_3d_ptr;

struct rmatrix_3d_builder
{
  private:
    std::size_t rows_;
    std::size_t columns_;
    std::size_t layers_;
    double dvalue_;

  public:
    LSS_API explicit rmatrix_3d_builder();

    LSS_API rmatrix_3d_builder &rows(std::size_t rows);

    LSS_API rmatrix_3d_builder &columns(std::size_t columns);

    LSS_API rmatrix_3d_builder &layers(std::size_t layers);

    LSS_API rmatrix_3d_builder &default_value(double value);

    LSS_API rmatrix_3d_ptr build();
};

struct cmatrix_3d_builder
{
  private:
    std::size_t rows_;
    std::size_t columns_;
    std::size_t layers_;
    double dvalue_;

  public:
    LSS_API explicit cmatrix_3d_builder();

    LSS_API cmatrix_3d_builder &rows(std::size_t rows);

    LSS_API cmatrix_3d_builder &columns(std::size_t columns);

    LSS_API cmatrix_3d_builder &layers(std::size_t layers);

    LSS_API cmatrix_3d_builder &default_value(double value);

    LSS_API cmatrix_3d_ptr build();
};

} // namespace lss

#endif ///_CONTAINER_2D_HPP_
