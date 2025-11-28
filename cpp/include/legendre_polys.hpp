#pragma once

#include <iostream>
#include <string>
#include <vector>

namespace flux {
class LegendrePolys {
    static void raise_error(const std::string &msg) {
        std::cerr << msg;
        exit(1);
    }

public:
    static double eval(std::size_t n, double x) {
        switch (n) {
        case 0: return f0(x);
        case 1: return f1(x);
        case 2: return f2(x);
        case 3: return f3(x);
        case 4: return f4(x);
        case 5: return f5(x);
        case 6: return f6(x);
        default:
            raise_error("LegendrePolys: out of range " + std::to_string(n));
            return 0;
        }
    }

    [[maybe_unused]] constexpr static double f0(double x) { return 1; };

    [[maybe_unused]] constexpr static double f1(double x) { return x; };

    [[maybe_unused]] constexpr static double f2(double x) {
        return (3 * x * x - 1) / 2.0;
    };

    [[maybe_unused]] constexpr static double f3(double x) {
        return (5 * x * x - 3) * x / 2.0;
    };

    [[maybe_unused]] constexpr static double f4(double x) {
        return ((35 * x * x - 30) * x * x + 3) / 8.0;
    };

    [[maybe_unused]] constexpr static double f5(double x) {
        return ((63 * x * x - 70) * x * x + 15) * x / 8.0;
    };

    [[maybe_unused]] constexpr static double f6(double x) {
        return (((231 * x * x - 315) * x * x + 105) * x * x - 5) / 16.0;
    };
};

class LegendrePolysDx {
    static void raise_error(const std::string &msg) {
        std::cerr << msg;
        exit(1);
    }

public:
    static double eval(std::size_t n, double x) {
        switch (n) {
        case 0: return f0(x);
        case 1: return f1(x);
        case 2: return f2(x);
        case 3: return f3(x);
        case 4: return f4(x);
        case 5: return f5(x);
        case 6: return f6(x);
        default:
            raise_error("LegendrePolysDx: out of range " + std::to_string(n));
            return 0;
        }
    }

    [[maybe_unused]] constexpr static double f0(double x) { return 0; };

    [[maybe_unused]] constexpr static double f1(double x) { return 1; };

    [[maybe_unused]] constexpr static double f2(double x) { return 3 * x; };

    [[maybe_unused]] constexpr static double f3(double x) {
        return (15 * x * x - 3) / 2.0;
    };

    [[maybe_unused]] constexpr static double f4(double x) {
        return (140 * x * x - 60) * x / 8.0;
    };

    [[maybe_unused]] constexpr static double f5(double x) {
        return ((315 * x * x - 210) * x * x + 15) / 8.0;
    };

    [[maybe_unused]] constexpr static double f6(double x) {
        return ((1386 * x * x - 1260) * x * x + 210) * x / 16.0;
    };
};

template <typename Poly>
double evals(const std::vector<double> &vec, double x, size_t id_start,
             size_t id_len) {
    double result = 0;
    for (size_t i = 0; i < id_len; i++) {
        result += (vec[id_start + i] * Poly::eval(i, x));
    }

    return result;
}
}  // namespace flux
