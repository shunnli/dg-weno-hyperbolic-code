#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gaussquad {

template <typename Func>
concept NumericalFunction = requires(Func f, double x, double w) {
    { f(x) * w } -> std::same_as<double>;
    { f(x) + f(x) } -> std::same_as<double>;
};

class Quad {
public:
    explicit Quad(const std::pair<std::vector<double>, std::vector<double>>
                      &points_and_weights)
        : m_points(points_and_weights.first),
          m_weights(points_and_weights.second),
          m_len(points_and_weights.first.size()) {
        if (m_points.size() != m_weights.size()) {
            throw std::runtime_error("points.size() != weights.size()");
        }
    }

    template <size_t N>
    explicit Quad(const std::pair<std::array<double, N>, std::array<double, N>>
                      &points_and_weights)
        : m_points(points_and_weights.first.begin(),
                   points_and_weights.first.end()),
          m_weights(points_and_weights.second.begin(),
                    points_and_weights.second.end()),
          m_len(N) {}

    const std::vector<double> &points() const { return m_points; }

    const std::vector<double> &weights() const { return m_weights; }

    template <NumericalFunction Func, NumericalFunction DFunc>
    Quad &transform(Func f, DFunc df) {
        std::transform(m_points.begin(), m_points.end(), m_points.begin(), f);
        std::transform(m_points.begin(), m_points.end(), m_weights.begin(),
                       m_weights.begin(),
                       [&df](auto x, auto w) { return df(x) * w; });
        return *this;
    }

    template <NumericalFunction Func>
    auto integrate(Func f) {
        return std::inner_product(
            m_points.cbegin(), m_points.cend(), m_weights.cbegin(), double{},
            std::plus<>(),
            [&f](double x, double w) -> double { return f(x) * w; });
    }

private:
    std::vector<double> m_points;
    std::vector<double> m_weights;
    std::size_t m_len;
};

}  // namespace gaussquad
