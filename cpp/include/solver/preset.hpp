#pragma once

#include <vector>

#include "requires.h"

namespace flux {

struct Mesh1d {
    double dx;
};

struct Vec {
    std::vector<double> data;

    explicit Vec(std::vector<double> d) : data(std::move(d)) {}

    Vec(const Vec &rhs) = default;

    Vec &operator=(const Vec &rhs) = default;

    Vec(Vec &&rhs) noexcept = default;

    Vec &operator=(Vec &&rhs) noexcept = default;

    ~Vec() = default;

    Vec operator+(const Vec &rhs) const {
        Vec result(data);
        for (size_t i = 0; i < data.size(); ++i) {
            result.data[i] += rhs.data[i];
        }
        return result;
    }

    Vec operator-(const Vec &rhs) const {
        Vec result(data);
        for (size_t i = 0; i < data.size(); ++i) {
            result.data[i] -= rhs.data[i];
        }
        return result;
    }

    friend Vec operator*(double scalar, const Vec &vec) {
        Vec result(vec);
        for (auto &v : result.data) { v *= scalar; }
        return result;
    }
};

static_assert(VarRequirements<Vec>, "Vec does not satisfy VarRequirements!");
}  // namespace flux
