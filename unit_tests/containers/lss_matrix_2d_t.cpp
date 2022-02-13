#include "lss_matrix_2d_t.hpp"

void basic_matrix_2d()
{
    auto const rows = 3;
    auto const cols = 4;

    matrix_2d mat(rows, cols);

    LSS_ASSERT(mat.columns() == cols, "Columns must be of same size");
    LSS_ASSERT(mat.rows() == rows, "Rows must be of same size");

    // 0.layer:
    mat(0, 0, 0);
    mat(0, 1, 1);
    mat(0, 2, 2);
    mat(0, 3, 3);
    //
    mat(1, 0, 4);
    mat(1, 1, 5);
    mat(1, 2, 6);
    mat(1, 3, 7);
    //
    mat(2, 0, 8);
    mat(2, 1, 9);
    mat(2, 2, 10);
    mat(2, 3, 11);

    print_matrix_2d(mat, "mat");

    auto data = mat.data();
    print_array(data, "data");
    print_matrix_2d(mat, "mat after data transfer");
    std::valarray<double> vals = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        LSS_ASSERT(data[i] == vals[i], "value at " << i << "must be identical");
    }

    auto shifted_data = vals.cshift(5);
    print_array(shifted_data, "shifted by 5 steps");
    mat.from_data(std::move(shifted_data));
    print_array(shifted_data, "shifted_data after populating matrix");
    LSS_ASSERT(shifted_data.size() == 0, "must be empty after moving data to matrix");
    print_matrix_2d(mat, "mat from data that were circularly shifted by 5 steps");
}

void slice_matrix_2d()
{
    auto const rows = 3;
    auto const cols = 4;

    matrix_2d mat(rows, cols);

    LSS_ASSERT(mat.columns() == cols, "Columns must be of same size");
    LSS_ASSERT(mat.rows() == rows, "Rows must be of same size");

    // 0.layer:
    mat(0, 0, 0);
    mat(0, 1, 1);
    mat(0, 2, 2);
    mat(0, 3, 3);
    //
    mat(1, 0, 4);
    mat(1, 1, 5);
    mat(1, 2, 6);
    mat(1, 3, 7);
    //
    mat(2, 0, 8);
    mat(2, 1, 9);
    mat(2, 2, 10);
    mat(2, 3, 11);

    print_matrix_2d(mat, "mat");

    // get the last row:
    const std::valarray<double> row_2 = mat.row(2);
    print_array(row_2, "row_2");
    std::valarray<double> rvals = {8, 9, 10, 11};

    for (std::size_t i = 0; i < row_2.size(); ++i)
    {
        LSS_ASSERT(row_2[i] == rvals[i], "value at " << i << "must be identical");
    }

    // get the 3rd column:
    const std::valarray<double> col_3 = mat.column(3);
    print_array(col_3, "col_3");
    std::valarray<double> cvals = {3, 7, 11};

    for (std::size_t i = 0; i < col_3.size(); ++i)
    {
        LSS_ASSERT(col_3[i] == cvals[i], "value at " << i << "must be identical");
    }

    // modify first row in place:
    mat.row(1) = {3.14, 6.28, 12.56, 25.12};
    print_matrix_2d(mat, "mat in place mod");
    rvals = {3.14, 6.28, 12.56, 25.12};

    for (std::size_t i = 0; i < row_2.size(); ++i)
    {
        LSS_ASSERT(mat(1, i) == rvals[i], "value at " << i << "must be identical");
    }

    // modify second row in place again:
    std::valarray<double> v = {0.14, 0.28, 1.56, 2.12};
    mat.row(2, std::move(v));
    print_matrix_2d(mat, "mat in place mod again");
    rvals = {0.14, 0.28, 1.56, 2.12};

    for (std::size_t i = 0; i < row_2.size(); ++i)
    {
        LSS_ASSERT(mat(2, i) == rvals[i], "value at " << i << "must be identical");
    }
}
