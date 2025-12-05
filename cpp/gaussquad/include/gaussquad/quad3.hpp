// clang-format off
/*
I(f) \approx I_n(f) = |\Delta| \sum_{i=1}^n w_i f(x_i, y_i, z_i)
where $(x,y,z)$ are barycentric coordinates.

|   i   |   w   |    x    |    y    |    z    |
| :---: | :---: | :-----: | :-----: | :-----: |
|   1   |  1.0  | 1.0 / 3 | 1.0 / 3 | 1.0 / 3 |

|   i   |    w    |    x    |    y    |    z    |       |
| :---: | :-----: | :-----: | :-----: | :-----: | :---: |
|   1   | 1.0 / 3 | 2.0 / 3 | 1.0 / 6 | 1.0 / 6 |  a1   |
|   2   | 1.0 / 3 | 2.0 / 3 | 1.0 / 6 | 1.0 / 6 |  a2   |
|   3   | 1.0 / 3 | 1.0 / 6 | 1.0 / 6 | 2.0 / 3 |  a3   |

|   i   |         w         |         x         |         y         |         z         |       |
| :---: | :---------------: | :---------------: | :---------------: | :---------------: | :---: |
|   1   |       0.225       |      1.0 / 3      |      1.0 / 3      |      1.0 / 3      |  a1   |
|   2   | 0.125939180544827 | 0.797426985353087 | 0.101286507323456 | 0.101286507323456 |  b1   |
|   3   | 0.125939180544827 | 0.101286507323456 | 0.797426985353087 | 0.101286507323456 |  b2   |
|   4   | 0.125939180544827 | 0.101286507323456 | 0.101286507323456 | 0.797426985353087 |  b3   |
|   5   | 0.132394152788506 | 0.059715871789770 | 0.470142064105115 | 0.470142064105115 |  c1   |
|   6   | 0.132394152788506 | 0.470142064105115 | 0.059715871789770 | 0.470142064105115 |  c2   |
|   7   | 0.132394152788506 | 0.470142064105115 | 0.470142064105115 | 0.059715871789770 |  c3   |

|   i   |         w         |         x         |         y         |         z         |       |
| :---: | :---------------: | :---------------: | :---------------: | :---------------: | :---: |
|   1   | 0.050844906370207 | 0.873821971016996 | 0.063089014491502 | 0.063089014491502 |  a1   |
|   2   | 0.050844906370207 | 0.063089014491502 | 0.873821971016996 | 0.063089014491502 |  a2   |
|   3   | 0.050844906370207 | 0.063089014491502 | 0.063089014491502 | 0.873821971016996 |  a3   |
|   4   | 0.116786275726379 | 0.501426509658179 | 0.249286745170910 | 0.249286745170911 |  b1   |
|   5   | 0.116786275726379 | 0.249286745170910 | 0.501426509658179 | 0.249286745170911 |  b2   |
|   6   | 0.116786275726379 | 0.249286745170910 | 0.249286745170911 | 0.501426509658179 |  b3   |
|   7   | 0.082851075618374 | 0.636502499121399 | 0.310352451033785 | 0.053145049844816 |  c1   |
|   8   | 0.082851075618374 | 0.636502499121399 | 0.053145049844816 | 0.310352451033785 |  c2   |
|   9   | 0.082851075618374 | 0.310352451033785 | 0.636502499121399 | 0.053145049844816 |  c3   |
|  10   | 0.082851075618374 | 0.310352451033785 | 0.053145049844816 | 0.636502499121399 |  c4   |
|  11   | 0.082851075618374 | 0.053145049844816 | 0.636502499121399 | 0.310352451033785 |  c5   |
|  12   | 0.082851075618374 | 0.053145049844816 | 0.310352451033785 | 0.636502499121399 |  c6   |
*/
// clang-format on

#pragma once

#include <array>
#include <concepts>
#include <stdexcept>
#include <vector>

namespace gaussquad {

template <typename Func>
concept NumericalFunction2D = requires(Func f, double x, double y, double w) {
    { f(x, y) * w } -> std::same_as<double>;
    { f(x, y) + f(x, y) } -> std::same_as<double>;
};

struct PointXY {
    double x;
    double y;
};

// Barycentric Coordinates
struct PointCoordinate {
    double p1;
    double p2;
    double p3;
};

struct Triangle {
    const double ax;
    const double ay;
    const double bx;
    const double by;
    const double cx;
    const double cy;

    PointXY trans_to_xy(double c1, double c2, double c3) const {
        double x = c1 * ax + c2 * bx + c3 * cx;
        double y = c1 * ay + c2 * by + c3 * cy;
        return PointXY{.x = x, .y = y};
    }

    PointCoordinate trans_to_coordinate(double x, double y) const {
        double p3 = ((ay - by) * x + (bx - ax) * y + ax * by - bx * ay)
                    / ((ay - by) * cx + (bx - ax) * cy + ax * by - bx * ay);

        double p2 = ((ay - cy) * x + (cx - ax) * y + ax * cy - cx * ay)
                    / ((ay - cy) * bx + (cx - ax) * by + ax * cy - cx * ay);

        double p1 = 1 - p2 - p3;

        return PointCoordinate{.p1 = p1, .p2 = p2, .p3 = p3};
    }

    double area() const { return triangle_signed_area(ax, ay, bx, by, cx, cy); }

    static double triangle_signed_area(double x1, double y1, double x2,
                                       double y2, double x3, double y3) {
        double px = x2 - x1;
        double py = y2 - y1;
        double qx = x3 - x1;
        double qy = y3 - y1;

        return (px * qy - py * qx) / 2;
    }
};

struct TriangleRule {
    static auto P1() -> std::pair<std::vector<double>, std::vector<double>> {
        std::vector<double> weights = {1.0};
        std::vector<double> points = {1.0 / 3, 1.0 / 3, 1.0 / 3};
        return std::make_pair(points, weights);
    }

    static auto P3() -> std::pair<std::vector<double>, std::vector<double>> {
        std::vector<double> weights = {1.0 / 3, 1.0 / 3, 1.0 / 3};
        std::vector<double> points = {
            2.0 / 3, 1.0 / 6, 1.0 / 6,  // p1
            1.0 / 6, 2.0 / 3, 1.0 / 6,  // p2
            1.0 / 6, 1.0 / 6, 2.0 / 3,  // p3
        };
        return std::make_pair(points, weights);
    }

    static auto P7() -> std::pair<std::vector<double>, std::vector<double>> {
        std::vector<double> weights = {
            0.225,  // 1
            0.125939180544827,
            0.125939180544827,
            0.125939180544827,  // 3
            0.132394152788506,
            0.132394152788506,
            0.132394152788506,  // 3
        };
        std::vector<double> points = {
            1.0 / 3, 1.0 / 3, 1.0 / 3,                                // p
                                                                      //
            0.797426985353087, 0.101286507323456, 0.101286507323456,  // q1
            0.101286507323456, 0.797426985353087, 0.101286507323456,  // q2
            0.101286507323456, 0.101286507323456, 0.797426985353087,  // q3
                                                                      //
            0.059715871789770, 0.470142064105115, 0.470142064105115,  // r1
            0.470142064105115, 0.059715871789770, 0.470142064105115,  // r2
            0.470142064105115, 0.470142064105115, 0.059715871789770,  // r3
        };
        return std::make_pair(points, weights);
    }

    static auto P12() -> std::pair<std::vector<double>, std::vector<double>> {
        std::vector<double> weights = {
            0.050844906370207, 0.050844906370207, 0.050844906370207,  // 3
            0.116786275726379, 0.116786275726379, 0.116786275726379,  // 3
            0.082851075618374, 0.082851075618374, 0.082851075618374,
            0.082851075618374, 0.082851075618374, 0.082851075618374,  // 6
        };
        std::vector<double> points = {
            0.873821971016996, 0.063089014491502, 0.063089014491502,  // p1
            0.063089014491502, 0.873821971016996, 0.063089014491502,  // p2
            0.063089014491502, 0.063089014491502, 0.873821971016996,  // p3
                                                                      //
            0.501426509658179, 0.249286745170910, 0.249286745170911,  // q1
            0.249286745170910, 0.501426509658179, 0.249286745170911,  // q2
            0.249286745170910, 0.249286745170911, 0.501426509658179,  // q3
                                                                      //
            0.636502499121399, 0.310352451033785, 0.053145049844816,  // r1
            0.636502499121399, 0.053145049844816, 0.310352451033785,  // r2
            0.310352451033785, 0.636502499121399, 0.053145049844816,  // r3
            0.310352451033785, 0.053145049844816, 0.636502499121399,  // r4
            0.053145049844816, 0.636502499121399, 0.310352451033785,  // r5
            0.053145049844816, 0.310352451033785, 0.636502499121399,  // r6
        };
        return std::make_pair(points, weights);
    }
};

class Quad3 {
public:
    explicit Quad3(const std::pair<std::vector<double>, std::vector<double>>
                       &points_and_weights)
        : m_points(points_and_weights.first),
          m_weights(points_and_weights.second),
          m_len(points_and_weights.second.size()) {
        if (m_points.size() != 3 * m_weights.size()) {
            throw std::runtime_error("points.size() != 3*weights.size()");
        }
    }

    template <size_t N>
    explicit Quad3(
        const std::pair<std::array<double, N>, std::array<double, 3 * N>>
            &points_and_weights)
        : m_points(points_and_weights.first.begin(),
                   points_and_weights.first.end()),
          m_weights(points_and_weights.second.begin(),
                    points_and_weights.second.end()),
          m_len(N) {}

    const std::vector<double> &points() const { return m_points; }

    const std::vector<double> &weights() const { return m_weights; }

    template <NumericalFunction2D Func>
    double integrate(Func f, const Triangle &the_triangle) const {
        double result = 0;
        for (std::size_t i = 0; i < m_len; ++i) {
            auto [x, y] = the_triangle.trans_to_xy(
                m_points[3 * i], m_points[3 * i + 1], m_points[3 * i + 2]);
            result += m_weights[i] * f(x, y);
        }
        return the_triangle.area() * result;
    }

private:
    std::vector<double> m_points;
    std::vector<double> m_weights;
    std::size_t m_len{};
};

}  // namespace gaussquad
