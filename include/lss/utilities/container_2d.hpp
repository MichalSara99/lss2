#if !defined(_CONTAINER_2D_HPP_)
#define _CONTAINER_2D_HPP_

#include <functional>

#include "../../../common/lss_macros.hpp"
#include "../../../common/lss_utility.hpp"
#include "../../../containers/lss_container_2d.hpp"

namespace lss
{

using rmatrix_2d = lss_containers::rmatrix_2d;
using rmatrix_2d_ptr = lss_containers::rmatrix_2d_ptr;
using cmatrix_2d = lss_containers::cmatrix_2d;
using cmatrix_2d_ptr = lss_containers::cmatrix_2d_ptr;

struct rmatrix_2d_builder
{
  private:
    std::size_t rows_;
    std::size_t columns_;
    double dvalue_;

  public:
    LSS_API explicit rmatrix_2d_builder();

    LSS_API rmatrix_2d_builder &rows(std::size_t rows);

    LSS_API rmatrix_2d_builder &columns(std::size_t columns);

    LSS_API rmatrix_2d_builder &default_value(double value);

    LSS_API rmatrix_2d_ptr build();
};

struct cmatrix_2d_builder
{
  private:
    std::size_t rows_;
    std::size_t columns_;
    double dvalue_;

  public:
    LSS_API explicit cmatrix_2d_builder();

    LSS_API cmatrix_2d_builder &rows(std::size_t rows);

    LSS_API cmatrix_2d_builder &columns(std::size_t columns);

    LSS_API cmatrix_2d_builder &default_value(double value);

    LSS_API cmatrix_2d_ptr build();
};

} // namespace lss

#endif ///_CONTAINER_2D_HPP_
