#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#if __has_include(<format>)
#include <format>
#define USE_FORMAT
#else
#include <iomanip>
#endif

namespace flux {

enum class ErrorType { Linf = 0, L1, L2 };

inline double error(const std::vector<double> &u1,
                    const std::vector<double> &u2, double dx,
                    ErrorType error_type) {
    size_t len = u1.size();
    if (u2.size() != len) {
        std::cerr << "error: u1 and u2 have different length" << std::endl;
        exit(1);
    }

    double result = 0;
    for (size_t i = 0; i < len; ++i) {
        double tmp = std::abs(u1[i] - u2[i]);

        if (error_type == ErrorType::Linf) {
            if (tmp > result) result = tmp;
        }

        if (error_type == ErrorType::L1) { result += tmp * dx; }

        if (error_type == ErrorType::L2) { result += tmp * tmp * dx; }
    }

    if (error_type == ErrorType::L2) { result = std::sqrt(result); }

    return result;
}

inline std::vector<double> order(const std::vector<double> &error,
                                 const std::vector<size_t> &nlist) {
    auto len = error.size();
    if (nlist.size() != len) {
        std::cerr << "order: Length of nlist must be equal to length of errors"
                  << std::endl;
        exit(1);
    }
    if (len <= 1) {
        std::cerr << "order: n must be greater than 1" << std::endl;
        exit(1);
    }

    std::vector<double> result(len);
    result[0] = 0;
    for (size_t k = 1; k < len; ++k) {
        result[k] = -(log(error[k - 1] / error[k])
                      / log(static_cast<double>(nlist[k - 1])
                            / static_cast<double>(nlist[k])));
    }

    return result;
}

#ifdef USE_FORMAT

inline void print_error_table(
    std::ostream &out, const std::vector<size_t> &nlist,
    const std::vector<double> &error_l1, const std::vector<double> &error_l2,
    const std::vector<double> &error_linf, const std::vector<double> &order_l1,
    const std::vector<double> &order_l2, const std::vector<double> &order_linf,
    char delimiter) {
    auto header = std::format("{:^5} {:c} {:^12} {:c} {:^8} {:c} {:^12} {:c} "
                              "{:^8} {:c} {:^12} {:c} {:^8}\n",
                              "n", delimiter, "error_1", delimiter, "order",
                              delimiter, "error_2", delimiter, "order",
                              delimiter, "error_inf", delimiter, "order");
    out << header;

    for (size_t i = 0; i < nlist.size(); ++i) {
        auto row = std::format(
            "{:^5} {:c} {:^12.2e} {:c} {:^8.2f} {:c} {:^12.2e} {:c} "
            "{:^8.2f} {:c} {:^12.2e} {:c} {:^8.2f}\n",
            nlist[i], delimiter, error_l1[i], delimiter, order_l1[i], delimiter,
            error_l2[i], delimiter, order_l2[i], delimiter, error_linf[i],
            delimiter, order_linf[i]);
        out << row;
    }

    out << std::endl;
    return;
}

#else

inline void print_error_table(
    std::ostream &out, const std::vector<size_t> &nlist,
    const std::vector<double> &error_l1, const std::vector<double> &error_l2,
    const std::vector<double> &error_linf, const std::vector<double> &order_l1,
    const std::vector<double> &order_l2, const std::vector<double> &order_linf,
    char delimiter) {
    out << std::setw(5) << "n" << delimiter << std::setw(12) << "error_1"
        << delimiter << std::setw(8) << "order" << delimiter << std::setw(12)
        << "error_2" << delimiter << std::setw(8) << "order" << delimiter
        << std::setw(12) << "error_inf" << delimiter << std::setw(8) << "order"
        << "\n";

    for (size_t i = 0; i < nlist.size(); ++i) {
        out << std::setw(5) << nlist[i] << delimiter << std::scientific
            << std::setw(12) << std::setprecision(2) << error_l1[i] << delimiter
            << std::fixed << std::setw(8) << std::setprecision(2) << order_l1[i]
            << delimiter << std::scientific << std::setw(12)
            << std::setprecision(2) << error_l2[i] << delimiter << std::fixed
            << std::setw(8) << std::setprecision(2) << order_l2[i] << delimiter
            << std::scientific << std::setw(12) << std::setprecision(2)
            << error_linf[i] << delimiter << std::fixed << std::setw(8)
            << std::setprecision(2) << order_linf[i] << "\n";
    }

    out << std::endl;
    return;
}

#endif

inline void print_error_table_to_file(
    const std::string &file_name, const std::vector<size_t> &nlist,
    const std::vector<double> &error_l1, const std::vector<double> &error_l2,
    const std::vector<double> &error_linf, const std::vector<double> &order_l1,
    const std::vector<double> &order_l2, const std::vector<double> &order_linf,
    char delimiter) {
    if (file_name.empty()) return;

    std::fstream f(file_name, std::ios::out);

    if (f.fail()) {
        std::cerr << "print_error_table_to_file: fail to open file "
                  << file_name << std::endl;
        exit(1);
    }

    print_error_table(f, nlist, error_l1, error_l2, error_linf, order_l1,
                      order_l2, order_linf, delimiter);

    f.close();

    std::cout << "print error table to file " << file_name << '\n';
    return;
}
}  // namespace flux
