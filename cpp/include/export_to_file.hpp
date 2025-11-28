#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace flux {
inline void export_to_file(const std::string &file_name,
                           const std::vector<double> &x,
                           const std::vector<double> &y, char delimiter) {
    if (file_name.empty()) return;

    std::fstream f(file_name, std::ios::out);

    if (f.fail()) {
        std::cerr << "export_to_file: fail to open file " << file_name
                  << std::endl;
        exit(1);
    }

    size_t row1 = x.size();
    size_t row2 = y.size();
    size_t row = (row1 < row2) ? row1 : row2;
    for (std::size_t i = 0; i < row; i++) {
        f << x[i] << delimiter << y[i] << "\n";
    }

    f.close();

    std::cout << "export to file " << file_name << '\n';
    return;
}

inline void export_to_file(const std::string &file_name,
                           const std::vector<double> &x,
                           const std::vector<double> &y,
                           const std::vector<double> &z, char delimiter) {
    if (file_name.empty()) return;

    std::fstream f(file_name, std::ios::out);

    if (f.fail()) {
        std::cerr << "export_to_file: fail to open file " << file_name
                  << std::endl;
        exit(1);
    }

    size_t row1 = x.size();
    size_t row2 = y.size();
    size_t row3 = z.size();
    size_t row = (row1 < row2) ? row1 : row2;
    row = (row < row3) ? row : row3;
    for (std::size_t i = 0; i < row; i++) {
        f << x[i] << delimiter << y[i] << delimiter << z[i] << "\n";
    }

    f.close();
    std::cout << "export to file " << file_name << '\n';
    return;
}
}  // namespace flux
