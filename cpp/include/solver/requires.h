#pragma once

template <typename T>
concept VarRequirements = requires(T a, T b, double d) {
    T(a + b);
    T(a - b);
    T(d * a);
    a = b;
    T(a);
};
