#pragma once

#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

// NOLINTBEGIN(readability-identifier-naming,google-explicit-constructor,hicpp-explicit-conversions)

namespace flux {

template <typename E>
class unexpected {
public:
    explicit unexpected(const E &err) : m_error(err) {}

    explicit unexpected(E &&err) : m_error(std::move(err)) {}

    const E &value() const & { return m_error; }

    E &value() & { return m_error; }

    E &&value() && { return std::move(m_error); }

private:
    E m_error;
};

template <typename T, typename E>
class expected {
    static_assert(!std::is_same_v<T, unexpected<E>>,
                  "T must not be unexpected<E>");

public:
    expected(const T &value) : m_data(value), m_valid(true) {}

    expected(T &&value) : m_data(std::move(value)), m_valid(true) {}

    expected(unexpected<E> err) : m_data(std::move(err)), m_valid(false) {}

    bool has_value() const noexcept { return m_valid; }

    explicit operator bool() const noexcept { return m_valid; }

    T &value() {
        if (!m_valid) throw std::runtime_error("Bad expected access: no value");
        return std::get<T>(m_data);
    }

    const T &value() const {
        if (!m_valid) throw std::runtime_error("Bad expected access: no value");
        return std::get<T>(m_data);
    }

    E &error() {
        if (m_valid) throw std::runtime_error("Bad expected access: no error");
        return std::get<unexpected<E>>(m_data).value();
    }

    const E &error() const {
        if (m_valid) throw std::runtime_error("Bad expected access: no error");
        return std::get<unexpected<E>>(m_data).value();
    }

private:
    std::variant<T, unexpected<E>> m_data;
    bool m_valid;
};

}  // namespace flux

// NOLINTEND(readability-identifier-naming,google-explicit-constructor,hicpp-explicit-conversions)
