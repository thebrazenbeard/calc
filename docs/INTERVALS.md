# Outward-rounded interval arithmetic

The `calc.intervals.Interval` module is a small, dependency-free implementation of **basic arithmetic on closed intervals with binary64 endpoints**. Its guarantee is intentionally narrow: for finite representable input endpoints, addition, subtraction, multiplication, and division (when the divisor excludes zero) enclose the corresponding exact real interval result.

```python
from calc.intervals import Interval

x = Interval.point(0.1)
y = Interval.point(0.2)
sum_bounds = x + y
assert sum_bounds.contains(0.1 + 0.2)
```

## Representation and guarantee

An interval is `[lower, upper]`, with ordered binary64 bounds. Integers are accepted only if they can be represented **exactly** by binary64; non-integer rational literals should be supplied as floating-point numbers with their already-rounded meaning made explicit.

For every operation, the implementation:

1. Converts finite binary64 bounds to **exact rational numbers** with `Fraction.from_float`.
2. Computes the exact rational endpoint extrema for the operation.
3. Converts those extrema back into binary64, stepping by `math.nextafter` toward negative or positive infinity *only when the nearest float would violate enclosure*.

This retains tight endpoints when the result is exactly representable and brackets rounding error otherwise. The guarantee covers the exact values of **input binary64 numbers**, not a user's intended decimal strings. For example, the machine number `0.1` is not the mathematical rational `1/10`.

| Operation | Exact lower extremum | Exact upper extremum |
| --- | --- | --- |
| `[a,b]+[c,d]` | `a+c` | `b+d` |
| `[a,b]-[c,d]` | `a-d` | `b-c` |
| `[a,b]*[c,d]` | Minimum of four endpoint products | Maximum of four products |
| `[a,b]/[c,d]` | Minimum of four endpoint quotients | Maximum of four quotients |

The quotient rule requires `0 ∉ [c,d]`.

## Domain and numerical boundaries

- If a real operation exceeds finite binary64 range, the enclosure can use `-inf` or `+inf` on one side. For example, `[1e308,1e308]*[2,2]` produces `[max_float,+inf]`; it does *not* silently return a falsely finite value.
- Arithmetic on an interval that already has an infinite endpoint is rejected rather than guessing how to handle indeterminate forms such as `0*inf`. Intersection, hull, negation, and containment remain defined for unbounded intervals.
- Subnormal underflow is enclosed via `[0,min_subnormal]` or `[-min_subnormal,-0]` as appropriate. Results equal to zero conventionally use a negative-zero lower bound and positive-zero upper bound.
- NaN and non-ordered endpoints are rejected. Division by an interval containing zero raises `ZeroDivisionError`.
- The implementation does **not** handle open intervals, complex intervals, dependency reduction, interval Newton methods, transcendental range bounding, arbitrary-precision rationals as *inputs*, or general validated ODE/PDE solutions.

### Why exact rational intermediates?

Using the host's default round-to-nearest arithmetic followed by a blind `nextafter` step is commonly conservative on IEEE-754 systems, but exact rational endpoints make the enclosure decision auditable and avoid reliance on global floating-point rounding modes. This is relatively expensive; the module is designed for teaching and small validation tasks, not high-throughput interval simulation.

## Verification

`tests/test_intervals.py` compares output bounds against independent exact rational arithmetic for hand-picked extremal cases and 150 deterministic seeded random finite input pairs. It includes underflow, overflow, signed zero, invalid values, denominator-zero rejection, intersection, and hull.

Source references:

- [Python `math.nextafter` documentation](https://docs.python.org/3.11/library/math.html#math.nextafter)
- [Python `fractions.Fraction` exact float conversion](https://docs.python.org/3.11/library/fractions.html#fractions.Fraction.from_float)
- [Python floating-point representation and error](https://docs.python.org/3.11/tutorial/floatingpoint.html)

No third-party source code was copied.
