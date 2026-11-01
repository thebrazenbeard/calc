"""Small, deterministic symbolic calculus over real-valued expressions.

This is not a general computer algebra system. Expressions have a deliberately
small node vocabulary; user source is parsed with an AST allowlist, never run.
Real logarithms, division, and negative powers retain their usual domains.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Mapping


Scalar = int | float | Fraction
_MAX_EXPONENT = 64
_MAX_INTEGER_BITS = 4096


class Expr:
    """Immutable syntax tree node supporting arithmetic and formal derivatives."""

    def eval(self, values: Mapping[str, Scalar]) -> float:
        raise NotImplementedError

    def diff(self, variable: str) -> Expr:
        if not isinstance(variable, str) or not variable.isidentifier():
            raise ValueError("differentiation variable must be an identifier")
        return self._differentiate(variable).simplify()

    def _differentiate(self, variable: str) -> Expr:
        raise NotImplementedError

    def simplify(self) -> Expr:
        """Perform safe local reductions, not general algebraic equivalence."""
        return self

    def __add__(self, other: Expr | Scalar) -> Expr:
        return _add(self, _coerce(other))

    def __radd__(self, other: Expr | Scalar) -> Expr:
        return _add(_coerce(other), self)

    def __sub__(self, other: Expr | Scalar) -> Expr:
        return _add(self, -_coerce(other))

    def __rsub__(self, other: Expr | Scalar) -> Expr:
        return _add(_coerce(other), -self)

    def __neg__(self) -> Expr:
        return _mul(Const(-1), self)

    def __pos__(self) -> Expr:
        return self

    def __mul__(self, other: Expr | Scalar) -> Expr:
        return _mul(self, _coerce(other))

    def __rmul__(self, other: Expr | Scalar) -> Expr:
        return _mul(_coerce(other), self)

    def __truediv__(self, other: Expr | Scalar) -> Expr:
        return _div(self, _coerce(other))

    def __rtruediv__(self, other: Expr | Scalar) -> Expr:
        return _div(_coerce(other), self)

    def __pow__(self, exponent: int) -> Expr:
        if isinstance(exponent, bool):
            raise TypeError("boolean exponent is not an integer power")
        if type(exponent) is not int:
            if isinstance(exponent, Const) and type(exponent.value) is int:
                exponent = exponent.value
            else:
                raise ValueError("only integer constant powers are supported")
        if abs(exponent) > _MAX_EXPONENT:
            raise ValueError("exponent exceeds supported resource bound")
        return _pow(self, exponent)


def _validated_scalar(value: Scalar) -> Scalar:
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction)):
        raise TypeError("only real int, float, or Fraction scalars are supported")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("non-finite numeric value")
    if isinstance(value, int) and value.bit_length() > _MAX_INTEGER_BITS:
        raise ValueError("integer exceeds supported resource bound")
    if isinstance(value, Fraction) and max(value.numerator.bit_length(), value.denominator.bit_length()) > _MAX_INTEGER_BITS:
        raise ValueError("rational exceeds supported resource bound")
    return value


def _coerce(value: Expr | Scalar) -> Expr:
    return value if isinstance(value, Expr) else Const(_validated_scalar(value))


def _float(value: Scalar) -> float:
    _validated_scalar(value)
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError("value too large for floating evaluation") from exc
    if not math.isfinite(result):
        raise ValueError("floating evaluation is non-finite")
    return result


def _finite(value: float) -> float:
    if not math.isfinite(value):
        raise ValueError("non-finite expression result")
    return value


@dataclass(frozen=True)
class Const(Expr):
    """Real constant. Integers and fractions remain exact during folding."""

    value: Scalar

    def __post_init__(self) -> None:
        _validated_scalar(self.value)

    def eval(self, values: Mapping[str, Scalar]) -> float:
        return _float(self.value)

    def _differentiate(self, variable: str) -> Expr:
        return Const(0)


@dataclass(frozen=True)
class Var(Expr):
    """Named independent variable."""

    name: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.isidentifier():
            raise ValueError("variable name must be a valid identifier")

    def eval(self, values: Mapping[str, Scalar]) -> float:
        return _float(values[self.name])

    def _differentiate(self, variable: str) -> Expr:
        return Const(1 if self.name == variable else 0)


@dataclass(frozen=True)
class Add(Expr):
    left: Expr
    right: Expr

    def eval(self, values: Mapping[str, Scalar]) -> float:
        return _finite(self.left.eval(values) + self.right.eval(values))

    def _differentiate(self, variable: str) -> Expr:
        return self.left.diff(variable) + self.right.diff(variable)

    def simplify(self) -> Expr:
        return _add(self.left.simplify(), self.right.simplify())


@dataclass(frozen=True)
class Mul(Expr):
    left: Expr
    right: Expr

    def eval(self, values: Mapping[str, Scalar]) -> float:
        return _finite(self.left.eval(values) * self.right.eval(values))

    def _differentiate(self, variable: str) -> Expr:
        return self.left.diff(variable) * self.right + self.left * self.right.diff(variable)

    def simplify(self) -> Expr:
        return _mul(self.left.simplify(), self.right.simplify())


@dataclass(frozen=True)
class Div(Expr):
    numerator: Expr
    denominator: Expr

    def eval(self, values: Mapping[str, Scalar]) -> float:
        return _finite(self.numerator.eval(values) / self.denominator.eval(values))

    def _differentiate(self, variable: str) -> Expr:
        a, b = self.numerator, self.denominator
        return (a.diff(variable) * b - a * b.diff(variable)) / (b**2)

    def simplify(self) -> Expr:
        return _div(self.numerator.simplify(), self.denominator.simplify())


@dataclass(frozen=True)
class Pow(Expr):
    base: Expr
    exponent: int

    def __post_init__(self) -> None:
        if type(self.exponent) is not int or abs(self.exponent) > _MAX_EXPONENT:
            raise ValueError("power must be a bounded integer exponent")

    def eval(self, values: Mapping[str, Scalar]) -> float:
        return _finite(self.base.eval(values) ** self.exponent)

    def _differentiate(self, variable: str) -> Expr:
        if self.exponent == 0:
            return Const(0)
        return Const(self.exponent) * self.base ** (self.exponent - 1) * self.base.diff(variable)

    def simplify(self) -> Expr:
        return _pow(self.base.simplify(), self.exponent)


@dataclass(frozen=True)
class Function(Expr):
    name: str
    argument: Expr

    def __post_init__(self) -> None:
        if self.name not in ("sin", "cos", "exp", "log"):
            raise ValueError("unsupported symbolic function")

    def eval(self, values: Mapping[str, Scalar]) -> float:
        value = self.argument.eval(values)
        if self.name == "log" and value <= 0:
            raise ValueError("real logarithm requires a positive argument")
        method = {"sin": math.sin, "cos": math.cos, "exp": math.exp, "log": math.log}[self.name]
        return _finite(method(value))

    def _differentiate(self, variable: str) -> Expr:
        inner = self.argument.diff(variable)
        if self.name == "sin":
            return cos(self.argument) * inner
        if self.name == "cos":
            return -sin(self.argument) * inner
        if self.name == "exp":
            return exp(self.argument) * inner
        if self.name == "log":
            return inner / self.argument
        raise ValueError("unsupported function")

    def simplify(self) -> Expr:
        arg = self.argument.simplify()
        if isinstance(arg, Const):
            return Const(Function(self.name, arg).eval({}))
        return Function(self.name, arg)


def _add(a: Expr, b: Expr) -> Expr:
    if isinstance(a, Const) and isinstance(b, Const):
        return Const(a.value + b.value)
    if a == Const(0):
        return b
    if b == Const(0):
        return a
    return Add(a, b)


def _mul(a: Expr, b: Expr) -> Expr:
    if isinstance(a, Const) and isinstance(b, Const):
        return Const(a.value * b.value)
    if a == Const(1):
        return b
    if b == Const(1):
        return a
    # Do not erase 0*f(x): f can be undefined at the evaluation point.
    return Mul(a, b)


def _div(a: Expr, b: Expr) -> Expr:
    if isinstance(a, Const) and isinstance(b, Const):
        if isinstance(a.value, (int, Fraction)) and isinstance(b.value, (int, Fraction)):
            return Const(Fraction(a.value) / Fraction(b.value))
        return Const(a.value / b.value)
    if b == Const(1):
        return a
    return Div(a, b)


def _pow(base: Expr, exponent: int) -> Expr:
    if exponent == 0:
        # The base may itself be undefined (for example (1/x)**0 at x=0).
        return Const(1) if isinstance(base, (Const, Var)) else Pow(base, 0)
    if exponent == 1:
        return base
    if isinstance(base, Const):
        try:
            if exponent < 0 and isinstance(base.value, (int, Fraction)):
                value = Fraction(base.value)**exponent
            else:
                value = base.value**exponent
        except OverflowError:
            return Pow(base, exponent)
        if isinstance(value, float) and not math.isfinite(value):
            return Pow(base, exponent)
        return Const(value)
    return Pow(base, exponent)


def sin(value: Expr | Scalar) -> Expr:
    return Function("sin", _coerce(value)).simplify()


def cos(value: Expr | Scalar) -> Expr:
    return Function("cos", _coerce(value)).simplify()


def exp(value: Expr | Scalar) -> Expr:
    return Function("exp", _coerce(value)).simplify()


def log(value: Expr | Scalar) -> Expr:
    return Function("log", _coerce(value)).simplify()


_FUNCTIONS = {"sin": sin, "cos": cos, "exp": exp, "log": log}


def parse_expression(source: str, *, max_nodes: int = 256, max_depth: int = 40) -> Expr:
    """Parse a bounded arithmetic expression, with no Python execution.

    Allowed: names, real numeric literals, + - * /, integer powers, unary +/-,
    and one-argument sin/cos/exp/log calls. No keywords, attributes,
    subscripts, comprehensions, lambdas, strings, assignments, or imports.
    """
    if not isinstance(source, str) or not source.strip() or len(source) > 2048:
        raise ValueError("expression must be a nonempty string of <= 2048 characters")
    if type(max_nodes) is not int or type(max_depth) is not int or max_nodes < 1 or max_depth < 1:
        raise ValueError("resource limits must be positive integers")
    try:
        tree = ast.parse(source, mode="eval")
    except SyntaxError as exc:
        raise ValueError("invalid expression syntax") from exc
    if sum(1 for _ in ast.walk(tree)) > max_nodes:
        raise ValueError("expression exceeds node limit")

    def read(node: ast.AST, depth: int) -> Expr:
        if depth > max_depth:
            raise ValueError("expression exceeds nesting limit")
        if isinstance(node, ast.Expression):
            return read(node.body, depth + 1)
        if isinstance(node, ast.Constant):
            if type(node.value) not in (int, float):
                raise ValueError("only real numeric literals are allowed")
            return Const(node.value)
        if isinstance(node, ast.Name):
            return Var(node.id)
        if isinstance(node, ast.UnaryOp):
            operand = read(node.operand, depth + 1)
            if isinstance(node.op, ast.UAdd):
                return operand
            if isinstance(node.op, ast.USub):
                return -operand
            raise ValueError("unsupported unary operator")
        if isinstance(node, ast.BinOp):
            left, right = read(node.left, depth + 1), read(node.right, depth + 1)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                if not isinstance(right, Const) or type(right.value) is not int:
                    raise ValueError("powers must be integer constants")
                return left**right.value
            raise ValueError("unsupported binary operator")
        if isinstance(node, ast.Call):
            if (not isinstance(node.func, ast.Name) or node.func.id not in _FUNCTIONS
                    or len(node.args) != 1 or node.keywords):
                raise ValueError("only one-argument sin/cos/exp/log are allowed")
            return _FUNCTIONS[node.func.id](read(node.args[0], depth + 1))
        raise ValueError(f"unsupported syntax: {node.__class__.__name__}")

    return read(tree, 0).simplify()


__all__ = ["Expr", "Const", "Var", "sin", "cos", "exp", "log", "parse_expression"]
