"""A bounded, pure *untyped lambda calculus* interpreter.

This module models term rewriting in computation theory. It has no connection
to the derivative of a real-valued function. It never executes Python code.

Terms: Var(name), Abs(parameter, body), App(function, argument).
Evaluation: normal-order (leftmost/outermost) beta reduction, including inside
lambda bodies, with alpha-renaming when substitution could capture variables.
Untyped terms may not have a normal form, hence mandatory step/size budgets.
"""

from __future__ import annotations

from dataclasses import dataclass
import re


_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z_0-9]*\Z")
_TOKEN = re.compile(r"\s*(\\|λ|\.|\(|\)|[A-Za-z_][A-Za-z_0-9]*)")
_SOURCE_LIMIT = 4096
_TOKEN_LIMIT = 1024
_PARSE_DEPTH_LIMIT = 128


class LambdaSyntaxError(ValueError):
    """A malformed or over-budget lambda expression."""


class ReductionLimitExceeded(RuntimeError):
    """Normal-order normalization did not finish within its resource limits."""


@dataclass(frozen=True, slots=True)
class Var:
    name: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not _IDENTIFIER.fullmatch(self.name):
            raise ValueError("variable must be a valid identifier")


@dataclass(frozen=True, slots=True)
class Abs:
    parameter: str
    body: Term

    def __post_init__(self) -> None:
        if not isinstance(self.parameter, str) or not _IDENTIFIER.fullmatch(self.parameter):
            raise ValueError("abstraction parameter must be a valid identifier")
        _validate_term(self.body)


@dataclass(frozen=True, slots=True)
class App:
    function: Term
    argument: Term

    def __post_init__(self) -> None:
        _validate_term(self.function)
        _validate_term(self.argument)


Term = Var | Abs | App


def _validate_term(term: object) -> None:
    if not isinstance(term, (Var, Abs, App)):
        raise TypeError("expected a lambda-calculus Var, Abs, or App")


def parse(source: str) -> Term:
    """Parse variables, applications, and ``\\x.term`` / ``λx.term``.

    Application associates left; an unparenthesized lambda abstraction extends
    rightward. For example ``(\\x.x) y`` is an application, while ``\\x.x y``
    is an abstraction whose body applies x to y. Only the specified grammar is
    accepted, not Python source or arbitrary function calls.
    """
    if not isinstance(source, str) or len(source) > _SOURCE_LIMIT:
        raise LambdaSyntaxError("source must be a string of at most 4096 characters")
    tokens: list[str] = []
    position = 0
    while position < len(source):
        if source[position:].isspace():
            break
        match = _TOKEN.match(source, position)
        if match is None:
            raise LambdaSyntaxError(f"unexpected character at offset {position}")
        tokens.append(match.group(1))
        if len(tokens) > _TOKEN_LIMIT:
            raise LambdaSyntaxError("too many tokens")
        position = match.end()
    if not tokens:
        raise LambdaSyntaxError("empty lambda term")

    cursor = 0

    def peek() -> str | None:
        return tokens[cursor] if cursor < len(tokens) else None

    def take() -> str:
        nonlocal cursor
        if cursor >= len(tokens):
            raise LambdaSyntaxError("unexpected end of lambda term")
        token = tokens[cursor]
        cursor += 1
        return token

    def require(symbol: str) -> None:
        if take() != symbol:
            raise LambdaSyntaxError(f"expected {symbol!r}")

    def parse_term(depth: int) -> Term:
        if depth > _PARSE_DEPTH_LIMIT:
            raise LambdaSyntaxError("lambda term is nested too deeply")
        if peek() in ("\\", "λ"):
            take()
            name = take()
            if not _IDENTIFIER.fullmatch(name):
                raise LambdaSyntaxError("expected a parameter name after lambda")
            require(".")
            return Abs(name, parse_term(depth + 1))

        result = parse_atom(depth + 1)
        # A bare lambda begins a whole term, so use parentheses when it is an
        # argument: f (\\x.x), never implicit Python or function-call syntax.
        while peek() == "(" or (peek() is not None and _IDENTIFIER.fullmatch(peek())):
            result = App(result, parse_atom(depth + 1))
        return result

    def parse_atom(depth: int) -> Term:
        if depth > _PARSE_DEPTH_LIMIT:
            raise LambdaSyntaxError("lambda term is nested too deeply")
        token = peek()
        if token == "(":
            take()
            result = parse_term(depth + 1)
            require(")")
            return result
        if token is not None and _IDENTIFIER.fullmatch(token):
            return Var(take())
        raise LambdaSyntaxError(f"expected variable or '('; got {token!r}")

    result = parse_term(0)
    if cursor != len(tokens):
        raise LambdaSyntaxError(f"unexpected trailing token {tokens[cursor]!r}")
    return result


def free_vars(term: Term) -> frozenset[str]:
    """Return names not bound by any enclosing abstraction."""
    _validate_term(term)
    if isinstance(term, Var):
        return frozenset((term.name,))
    if isinstance(term, Abs):
        return free_vars(term.body) - {term.parameter}
    return free_vars(term.function) | free_vars(term.argument)


def _all_names(term: Term) -> frozenset[str]:
    if isinstance(term, Var):
        return frozenset((term.name,))
    if isinstance(term, Abs):
        return _all_names(term.body) | {term.parameter}
    return _all_names(term.function) | _all_names(term.argument)


def _fresh(base: str, unavailable: frozenset[str]) -> str:
    index = 0
    while f"{base}_{index}" in unavailable:
        index += 1
    return f"{base}_{index}"


def substitute(term: Term, variable: str, replacement: Term) -> Term:
    """Capture-avoiding substitution ``term[variable := replacement]``.

    A binder is alpha-renamed only if a free variable of ``replacement`` would
    otherwise become bound during substitution. Fresh names avoid *all* names
    in the relevant terms, not merely their free variables.
    """
    _validate_term(term)
    _validate_term(replacement)
    if not isinstance(variable, str) or not _IDENTIFIER.fullmatch(variable):
        raise ValueError("substitution variable must be a valid identifier")
    if variable not in free_vars(term):
        return term
    if isinstance(term, Var):
        return replacement if term.name == variable else term
    if isinstance(term, App):
        return App(substitute(term.function, variable, replacement),
                   substitute(term.argument, variable, replacement))
    if term.parameter == variable:
        return term
    if term.parameter not in free_vars(replacement):
        return Abs(term.parameter, substitute(term.body, variable, replacement))
    unavailable = _all_names(term) | _all_names(replacement) | {variable}
    fresh = _fresh(term.parameter, unavailable)
    renamed_body = substitute(term.body, term.parameter, Var(fresh))
    return Abs(fresh, substitute(renamed_body, variable, replacement))


def beta_step(term: Term) -> Term | None:
    """Perform one leftmost/outermost beta step; return None at normal form."""
    _validate_term(term)
    if isinstance(term, Var):
        return None
    if isinstance(term, Abs):
        reduced = beta_step(term.body)
        return None if reduced is None else Abs(term.parameter, reduced)
    if isinstance(term.function, Abs):
        return substitute(term.function.body, term.function.parameter, term.argument)
    function = beta_step(term.function)
    if function is not None:
        return App(function, term.argument)
    argument = beta_step(term.argument)
    return None if argument is None else App(term.function, argument)


def _count_nodes(
    term: Term, *, limit: int | None = None, depth_limit: int | None = None,
) -> int:
    """Count without recursion; reject deep manually constructed terms."""
    count = 0
    pending = [(term, 1)]
    while pending:
        node, depth = pending.pop()
        count += 1
        if limit is not None and count > limit:
            return count
        if depth_limit is not None and depth > depth_limit:
            raise ReductionLimitExceeded("term nesting exceeds supported depth")
        if isinstance(node, Abs):
            pending.append((node.body, depth + 1))
        elif isinstance(node, App):
            pending.extend(((node.function, depth + 1),
                            (node.argument, depth + 1)))
    return count


def _free_occurrences(term: Term, name: str) -> int:
    """Count replaceable variables, ignoring binders that shadow ``name``."""
    if isinstance(term, Var):
        return int(term.name == name)
    if isinstance(term, Abs):
        return 0 if term.parameter == name else _free_occurrences(term.body, name)
    return _free_occurrences(term.function, name) + _free_occurrences(term.argument, name)


def _next_reduction_size_change(term: Term) -> int | None:
    """Compute the exact node-count change of the next normal-order beta step.

    This calculation occurs before allocating the reduct. Alpha-renaming
    preserves node counts, and each substituted free occurrence replaces one
    variable node with one copy of the argument tree.
    """
    if isinstance(term, Var):
        return None
    if isinstance(term, Abs):
        return _next_reduction_size_change(term.body)
    if isinstance(term.function, Abs):
        argument_size = _count_nodes(term.argument)
        copies = _free_occurrences(term.function.body, term.function.parameter)
        return copies * (argument_size - 1) - 2 - argument_size
    first = _next_reduction_size_change(term.function)
    return first if first is not None else _next_reduction_size_change(term.argument)


def normalize(term: Term, *, max_steps: int = 1000, max_nodes: int = 20000) -> Term:
    """Reduce to beta-normal form, or raise when bounded resources are spent.

    Normal order can find a beta-normal form when one exists, but untyped terms
    may diverge. ``max_steps`` counts beta reductions and ``max_nodes`` bounds
    intermediate AST size; neither is a proof of termination.
    """
    _validate_term(term)
    if type(max_steps) is not int or max_steps < 0:
        raise ValueError("max_steps must be a nonnegative integer")
    if type(max_nodes) is not int or max_nodes < 1:
        raise ValueError("max_nodes must be a positive integer")
    node_count = _count_nodes(term, limit=max_nodes, depth_limit=256)
    if node_count > max_nodes:
        raise ReductionLimitExceeded("initial term exceeds node budget")
    try:
        for _ in range(max_steps):
            change = _next_reduction_size_change(term)
            if change is None:
                return term
            if node_count + change > max_nodes:
                raise ReductionLimitExceeded("next reduction would exceed node budget")
            reduced = beta_step(term)
            if reduced is None:
                raise AssertionError("normal-order step disagrees with size preflight")
            node_count += change
            term = reduced
        if _next_reduction_size_change(term) is not None:
            raise ReductionLimitExceeded("beta-reduction step budget exceeded")
    except RecursionError as exc:
        raise ReductionLimitExceeded("term nesting exceeds supported recursion depth") from exc
    return term


def _binding_key(term: Term, environment: tuple[str, ...]) -> tuple:
    if isinstance(term, Var):
        for distance, bound_name in enumerate(reversed(environment)):
            if term.name == bound_name:
                return ("bound", distance)
        return ("free", term.name)
    if isinstance(term, Abs):
        return ("abs", _binding_key(term.body, environment + (term.parameter,)))
    return ("app", _binding_key(term.function, environment),
            _binding_key(term.argument, environment))


def alpha_equivalent(left: Term, right: Term) -> bool:
    """Compare binder structure, ignoring names of bound variables."""
    _validate_term(left)
    _validate_term(right)
    return _binding_key(left, ()) == _binding_key(right, ())


def format_term(term: Term) -> str:
    """Format a term in a parseable, parenthesized human-readable notation."""
    _validate_term(term)
    if isinstance(term, Var):
        return term.name
    if isinstance(term, Abs):
        return f"λ{term.parameter}.{format_term(term.body)}"
    function = format_term(term.function)
    argument = format_term(term.argument)
    if isinstance(term.function, Abs):
        function = f"({function})"
    if isinstance(term.argument, (Abs, App)):
        argument = f"({argument})"
    return f"{function} {argument}"
