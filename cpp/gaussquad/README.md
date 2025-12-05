# gaussquad

## Gauss-Legendre & Gauss-Lobatto Quadrature

This library computes the nodes $\{x_i\}$ and weights $\{w_i\}$ for Gauss–Legendre and Gauss–Lobatto quadrature rules on the standard interval $[-1,1]$:

$$
I(f) = \int_{-1}^1 f(x) dx \approx I_n(f) = \sum_{i=1}^n w_i f(x_i).
$$

Function Declarations

```cpp
// Compile-time version
template <unsigned N>
consteval auto gausslegendre()
    -> std::pair<std::array<double, N>, std::array<double, N>>;

template <unsigned N>
consteval auto gausslobatto()
    -> std::pair<std::array<double, N>, std::array<double, N>>;

// Runtime version
inline auto gausslegendre(unsigned n)
    -> std::pair<std::vector<double>, std::vector<double>>;

inline auto gausslobatto(unsigned n)
    -> std::pair<std::vector<double>, std::vector<double>>;
```

Examples

```cpp
// Compile-time quadrature rules
auto [points, weights] = gausslegendre<3>();
auto [points, weights] = gausslobatto<4>();

// Runtime quadrature rules
auto [points, weights] = gausslegendre(5);
auto [points, weights] = gausslobatto(6);
```

References

- [Legendre-Gauss Quadrature Weights and Nodes](https://ww2.mathworks.cn/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes?s_tid=srchtitle_support_results_4_Gauss%20Lobatto)
- [Legendre-Gauss-Lobatto Nodes and Weights](https://ww2.mathworks.cn/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights?s_tid=srchtitle_support_results_3_Gauss%2520Lobatto)


## Quadrature Wrapper

The `Quad` class provides a convenient interface for evaluating integrals using the computed nodes $\{x_i\}$ and weights $\{w_i\}$.

$$
I(f) = \int_{-1}^1 f(x) dx \approx \sum_{i=1}^n w_i f(x_i).
$$

Example
```cpp
double result = Quad(gausslegendre(5)).integrate([](double x) { return x*x; });
std::cout << "NInt(x^2,{x,-1,1}) = " << result << "\n";
```


## Interval Transform

Quadrature on a general interval $\Omega = [a,b]$ can be obtained by a change of variables $y = h(x)$ mapping $[-1,1] \to \Omega$:

$$
\int_\Omega g(y)dy
= \int_{-1}^1 g(h(x)) h'(x) dx
\approx \sum_{i=1}^n w_i h'(x_i) g(h(x_i)).
$$

Example: Integration $g(y) = \sin(y)$ on $[0,\pi]$.

```cpp
// Mapping [-1,1] -> [0,pi]
auto pi_f  = [](double x) { return std::numbers::pi / 2 * (x + 1); };
auto pi_df = [](double)   { return std::numbers::pi / 2; };

double result = Quad(gausslegendre(5))
                    .transform(pi_f, pi_df)
                    .integrate([](double x) { return std::sin(x); });
std::cout << "NInt(sin(x),{x,0,pi}) = " << result << "\n";
```
