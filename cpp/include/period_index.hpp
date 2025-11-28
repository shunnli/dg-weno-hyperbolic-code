#pragma once

#include <cstddef>
namespace flux {
class PeriodIndex {
public:
    explicit PeriodIndex(size_t n, size_t i) : m_n(n), m_i(i) {}

    size_t c() const { return m_i; }

    size_t l() const { return (m_i + m_n - 1) % m_n; }

    size_t l(size_t step) const { return (m_i + m_n - step) % m_n; }

    size_t r() const { return (m_i + 1) % m_n; }

    size_t r(size_t step) const { return (m_i + step) % m_n; }

private:
    size_t m_n{0};
    size_t m_i{0};
};
}  // namespace flux
