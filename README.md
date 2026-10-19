# Adaptative quadratures implementation

Numerical Analysis assignment on numerical integration: simple Newton-Cotes quadratures,
composite quadratures (trapezoidal, Simpson and Gauss) and an adaptive Simpson quadrature,
applied to the integral of `f(x) = sin(exp(2x))` on `[0, 2]`.

The assignment report (in Catalan) is [`Quadratures_adaptatives_PolTibau.pdf`](Quadratures_adaptatives_PolTibau.pdf).

## Structure

| File | Contents |
| --- | --- |
| `Simple_quadratures.py` | `trobaPunts` and `trobaPesos`: equispaced nodes and Newton-Cotes weights obtained by solving the moment system. Also contains, commented out, the experiments with `exp(-x^2)` and `1/(1+16x^2)` |
| `Composite_quadratures.py` | `trapeziCompost`, `SimpsonCompost` and `GaussCompost`, plus the `calcula` function that applies weights and nodes |
| `Adaptative_quadratures.py` | Exercises 1, 2, 4 and 5: errors of the composite trapezoidal and Simpson rules, number of subintervals needed, and `SimpsonAdaptativa` |
| `quadratures_Gauss.py` | `quadraturaGauss(nIP)`: Gauss-Legendre nodes and weights from 1 to 15 points on `[-1, 1]`. Provided by the course |

## Running

```
pip install -r requirements.txt
python Adaptative_quadratures.py
```

The script prints the approximations and errors and shows the assignment plots. It can take
a while, because exercise 2 increases the number of subintervals one by one until the
tolerance is reached.
