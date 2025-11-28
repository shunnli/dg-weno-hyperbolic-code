#pragma once

#include <vector>

namespace flux {
inline std::vector<double> linespace(double start, double end, size_t num,
                                     bool end_point, double &ret_step) {
    std::vector<double> result;

    if (end_point) {
        if (num <= 1) return result;
        ret_step = (end - start) / static_cast<double>(num - 1);
    }
    else {
        if (num == 0) return result;
        ret_step = (end - start) / static_cast<double>(num);
    }

    for (size_t i = 0; i < num; ++i) {
        result.push_back(start + static_cast<double>(i) * ret_step);
    }

    return result;
}

inline std::vector<double> linespace_mid(double start, double end, size_t num,
                                         double &ret_step) {
    std::vector<double> result;

    if (num == 0) return result;
    ret_step = (end - start) / static_cast<double>(num);

    for (size_t i = 0; i < num; ++i) {
        result.push_back(start + (static_cast<double>(i) + 0.5) * ret_step);
    }

    return result;
}
}  // namespace flux
