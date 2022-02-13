#include "lss_matrix_3d_t.hpp"
#include "lss_matrix_2d_t.hpp"

void basic_matrix_3d()
{

    auto const rows = 3;
    auto const cols = 4;
    auto const lays = 2;

    matrix_3d mat(rows, cols, lays);

    LSS_ASSERT(mat.columns() == cols, "Columns must be of same size");
    LSS_ASSERT(mat.rows() == rows, "Rows must be of same size");
    LSS_ASSERT(mat.layers() == lays, "Layers must be of same size");

    // 0.layer:
    mat(0, 0, 0, 0);
    mat(0, 1, 0, 1);
    mat(0, 2, 0, 2);
    mat(0, 3, 0, 3);
    //
    mat(1, 0, 0, 4);
    mat(1, 1, 0, 5);
    mat(1, 2, 0, 6);
    mat(1, 3, 0, 7);
    //
    mat(2, 0, 0, 8);
    mat(2, 1, 0, 9);
    mat(2, 2, 0, 10);
    mat(2, 3, 0, 11);

    // 1.layer:
    mat(0, 0, 1, 12);
    mat(0, 1, 1, 13);
    mat(0, 2, 1, 14);
    mat(0, 3, 1, 15);
    //
    mat(1, 0, 1, 16);
    mat(1, 1, 1, 17);
    mat(1, 2, 1, 18);
    mat(1, 3, 1, 19);
    //
    mat(2, 0, 1, 20);
    mat(2, 1, 1, 21);
    mat(2, 2, 1, 22);
    mat(2, 3, 1, 23);

    print_matrix_3d(mat, "mat");

    auto const data = mat.data();
    print_array(data, "data");
    std::valarray<double> vals = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        LSS_ASSERT(data[i] == vals[i], "value at " << i << "must be identical");
    }

    auto shifted_data = vals.cshift(5);
    print_array(shifted_data, "shifted by 5 steps");
    mat.from_data(std::move(shifted_data));
    LSS_ASSERT(shifted_data.size() == 0, "must by empty after transfer of data");
    print_matrix_3d(mat, "mat from data that were circularly shifted by 5 steps");
}

void slice_matrix_3d()
{
    auto const rows = 3;
    auto const cols = 4;
    auto const lays = 2;

    matrix_3d mat(rows, cols, lays);

    LSS_ASSERT(mat.columns() == cols, "Columns must be of same size");
    LSS_ASSERT(mat.rows() == rows, "Rows must be of same size");
    LSS_ASSERT(mat.layers() == lays, "Layers must be of same size");

    // 0.layer:
    mat(0, 0, 0) = 0;
    mat(0, 1, 0) = 1;
    mat(0, 2, 0) = 2;
    mat(0, 3, 0) = 3;
    //
    mat(1, 0, 0) = 4;
    mat(1, 1, 0) = 5;
    mat(1, 2, 0) = 6;
    mat(1, 3, 0) = 7;
    //
    mat(2, 0, 0) = 8;
    mat(2, 1, 0) = 9;
    mat(2, 2, 0) = 10;
    mat(2, 3, 0) = 11;

    // 1.layer:
    mat(0, 0, 1) = 12;
    mat(0, 1, 1) = 13;
    mat(0, 2, 1) = 14;
    mat(0, 3, 1) = 15;
    //
    mat(1, 0, 1) = 16;
    mat(1, 1, 1) = 17;
    mat(1, 2, 1) = 18;
    mat(1, 3, 1) = 19;
    //
    mat(2, 0, 1) = 20;
    mat(2, 1, 1) = 21;
    mat(2, 2, 1) = 22;
    mat(2, 3, 1) = 23;

    print_matrix_3d(mat, "mat");

    auto const data = mat.data();
    print_array(data, "data");
    std::valarray<double> vals = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        LSS_ASSERT(data[i] == vals[i], "value at " << i << "must be identical");
    }

    // get the row at (2, 1):
    const std::valarray<double> row_21 = mat.row(2, 1);
    print_array(row_21, "row_21");
    std::valarray<double> rvals = {20, 21, 22, 23};

    for (std::size_t i = 0; i < row_21.size(); ++i)
    {
        LSS_ASSERT(row_21[i] == rvals[i], "value at " << i << "must be identical");
    }

    // get the column at (2, 1):
    const std::valarray<double> col_21 = mat.column(2, 1);
    print_array(col_21, "col_21");
    std::valarray<double> cvals = {14, 18, 22};

    for (std::size_t i = 0; i < col_21.size(); ++i)
    {
        LSS_ASSERT(col_21[i] == cvals[i], "value at " << i << "must be identical");
    }

    // get the layer at (2, 3):
    const std::valarray<double> lay_23 = mat.layer(2, 3);
    print_array(lay_23, "lay_23");
    std::valarray<double> lvals = {11, 23};

    for (std::size_t i = 0; i < lay_23.size(); ++i)
    {
        LSS_ASSERT(lay_23[i] == lvals[i], "value at " << i << "must be identical");
    }

    // modify row at (2, 1) in place:
    mat.row(2, 1) = {3.14, 6.28, 12.56, 25.12};
    print_matrix_3d(mat, "mat in place mod");
    rvals = {3.14, 6.28, 12.56, 25.12};

    for (std::size_t i = 0; i < rvals.size(); ++i)
    {
        LSS_ASSERT(mat(2, i, 1) == rvals[i], "value at " << i << "must be identical");
    }

    // modify second row  again:
    std::valarray<double> v({0.14, 0.28, 1.56, 2.12});
    mat.row(2, 0) = std::move(v);
    print_matrix_3d(mat, "mat mod again");
    rvals = {0.14, 0.28, 1.56, 2.12};

    for (std::size_t i = 0; i < rvals.size(); ++i)
    {
        LSS_ASSERT(mat(2, i, 0) == rvals[i], "value at " << i << "must be identical");
    }

    // modify column at (2, 1) in place:
    mat.column(2, 1) = {30, 60, 12};
    print_matrix_3d(mat, "mat in place mod");
    cvals = {30, 60, 12};

    for (std::size_t i = 0; i < cvals.size(); ++i)
    {
        LSS_ASSERT(mat(i, 2, 1) == cvals[i], "value at " << i << "must be identical");
    }

    // modify column at (2, 0) again:
    std::valarray<double> c({20, 60, 10});
    mat.column(2, 0) = std::move(c);
    print_matrix_3d(mat, "mat mod again");
    rvals = {20, 60, 10};

    for (std::size_t i = 0; i < rvals.size(); ++i)
    {
        LSS_ASSERT(mat(i, 2, 0) == rvals[i], "value at " << i << "must be identical");
    }

    // modify layer at (2, 1) in place:
    mat.layer(2, 1) = {300, 600};
    print_matrix_3d(mat, "mat in place mod");
    cvals = {300, 600};

    for (std::size_t i = 0; i < cvals.size(); ++i)
    {
        LSS_ASSERT(mat(2, 1, i) == cvals[i], "value at " << i << "must be identical");
    }

    // modify layer at (2, 0) again:
    std::valarray<double> l({200, 600});
    mat.layer(2, 0) = std::move(l);
    print_matrix_3d(mat, "mat mod again");
    rvals = {200, 600};

    for (std::size_t i = 0; i < rvals.size(); ++i)
    {
        LSS_ASSERT(mat(2, 0, i) == rvals[i], "value at " << i << "must be identical");
    }
}
