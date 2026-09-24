"""Numerically evaluate RGBeta ``report_betas`` expressions.

The Wolfram runner exports the reporting convention

    16*pi^2 dX/dln(mu) = beta_X^(1),

with gauge beta functions already converted from RGBeta's native
``d(g^2)/dln(mu)`` convention.  This module evaluates those exported
``InputForm`` strings directly on a canonical numerical UV state.

Supported RGBeta structures are the ones emitted by the current T3 runner:

    Bar[x]                      -> complex conjugation
    Trans[x]                    -> transpose
    Tr[a . b . ...]             -> matrix trace
    Matrix[a,b,...][indices]     -> matrix product, with displayed indices
                                    carrying no extra numerical operation

Ordinary arithmetic ``+ - * / ^`` is also supported.  Parsing is performed by
a small dedicated parser; Python ``eval`` is intentionally not used.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import re
from typing import Any, Mapping

import numpy as np

from Numerical.BetaVector import beta_values_to_derivative, validate_beta_keys
from Numerical.State import SharedT3UVState, T3UVState
from Numerical.StateVector import StateVectorLayout, pack_uv_state


Number = float | complex
Value = Number | np.ndarray


_TOKEN_RE = re.compile(
    r"""
    (?P<WS>\s+)
    |(?P<NUMBER>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)
    |(?P<NAME>[A-Za-z_$][A-Za-z0-9_$]*)
    |(?P<OP>[\+\-\*/\^\.,\(\)\[\]])
    """,
    re.VERBOSE,
)


@dataclass(frozen=True)
class _Token:
    kind: str
    text: str


def _tokenize(text: str) -> list[_Token]:
    tokens: list[_Token] = []
    position = 0

    while position < len(text):
        match = _TOKEN_RE.match(text, position)
        if match is None:
            snippet = text[position:position + 40]
            raise ValueError(
                "Unsupported RGBeta InputForm syntax near "
                f"{snippet!r}."
            )

        position = match.end()
        kind = match.lastgroup
        if kind == "WS":
            continue

        tokens.append(_Token(kind=kind or "", text=match.group()))

    tokens.append(_Token(kind="EOF", text=""))
    return tokens


def _is_array(value: Value) -> bool:
    return isinstance(value, np.ndarray)


def _ordinary_multiply(left: Value, right: Value) -> Value:
    """Mathematica Times as used by the current exports.

    Matrix-matrix products in RGBeta are emitted explicitly as ``Matrix[...]``
    or with ``.`` inside traces.  Therefore ``*`` here is scalar
    multiplication; multiplying two arrays with ``*`` is rejected rather than
    silently interpreted elementwise.
    """

    if _is_array(left) and _is_array(right):
        raise ValueError(
            "RGBeta expression used '*' between two matrices. "
            "The current evaluator requires Matrix[...] or '.' for "
            "matrix multiplication."
        )
    return left * right


def _matrix_multiply(left: Value, right: Value) -> Value:
    if not (_is_array(left) and _is_array(right)):
        raise ValueError(
            "RGBeta '.' product currently requires matrix operands."
        )
    return left @ right


class _Parser:
    def __init__(self, text: str, environment: Mapping[str, Any]) -> None:
        self.tokens = _tokenize(text)
        self.index = 0
        self.environment = environment

    @property
    def current(self) -> _Token:
        return self.tokens[self.index]

    def advance(self) -> _Token:
        token = self.current
        self.index += 1
        return token

    def accept(self, text: str) -> bool:
        if self.current.text == text:
            self.advance()
            return True
        return False

    def expect(self, text: str) -> None:
        if not self.accept(text):
            raise ValueError(
                f"Expected {text!r}, found {self.current.text!r}."
            )

    def parse(self) -> Value:
        value = self.parse_expression()
        if self.current.kind != "EOF":
            raise ValueError(
                f"Unexpected trailing token {self.current.text!r}."
            )
        return value

    def parse_expression(self) -> Value:
        value = self.parse_term()

        while self.current.text in {"+", "-"}:
            operator = self.advance().text
            rhs = self.parse_term()
            value = value + rhs if operator == "+" else value - rhs

        return value

    def parse_term(self) -> Value:
        value = self.parse_power()

        while self.current.text in {"*", "/", "."}:
            operator = self.advance().text
            rhs = self.parse_power()

            if operator == "*":
                value = _ordinary_multiply(value, rhs)
            elif operator == "/":
                if _is_array(rhs):
                    raise ValueError(
                        "Division by a matrix is not supported."
                    )
                value = value / rhs
            else:
                value = _matrix_multiply(value, rhs)

        return value

    def parse_power(self) -> Value:
        value = self.parse_unary()

        if self.accept("^"):
            exponent = self.parse_power()
            if _is_array(exponent):
                raise ValueError("Matrix exponent is not supported.")
            if _is_array(value):
                raise ValueError(
                    "Elementwise/matrix powers are not used by the current "
                    "RGBeta T3 export."
                )
            value = value ** exponent

        return value

    def parse_unary(self) -> Value:
        if self.accept("+"):
            return self.parse_unary()
        if self.accept("-"):
            return -self.parse_unary()
        return self.parse_primary()

    def parse_primary(self) -> Value:
        token = self.current

        if token.kind == "NUMBER":
            self.advance()
            value: Value = float(token.text)
            return value

        if token.text == "(":
            self.advance()
            value = self.parse_expression()
            self.expect(")")
            return value

        if token.kind == "NAME":
            name = self.advance().text

            if self.accept("["):
                args = self.parse_arguments_after_open_bracket()
                value = self.call_function(name, args)
            else:
                if name not in self.environment:
                    raise ValueError(
                        f"Unknown RGBeta symbol {name!r}."
                    )
                value = self.environment[name]

            # RGBeta prints free matrix indices as a second bracket:
            # Matrix[...][gen[$i], gen[$j]].  The numerical object returned by
            # Matrix[...] is already the full matrix, so consume and ignore the
            # displayed component indices.
            while self.accept("["):
                self.skip_index_arguments_after_open_bracket()

            return value

        raise ValueError(
            f"Unexpected RGBeta token {token.text!r}."
        )

    def parse_arguments_after_open_bracket(self) -> list[Value]:
        args: list[Value] = []

        if self.accept("]"):
            return args

        while True:
            args.append(self.parse_expression())
            if self.accept("]"):
                break
            self.expect(",")

        return args

    def skip_index_arguments_after_open_bracket(self) -> None:
        depth = 1

        while depth:
            token = self.advance()
            if token.kind == "EOF":
                raise ValueError("Unterminated RGBeta index bracket.")
            if token.text == "[":
                depth += 1
            elif token.text == "]":
                depth -= 1

    def call_function(self, name: str, args: list[Value]) -> Value:
        if name == "Bar":
            self._require_arity(name, args, 1)
            return np.conjugate(args[0])

        if name == "Trans":
            self._require_arity(name, args, 1)
            value = args[0]
            if not _is_array(value):
                return value
            return np.transpose(value)

        if name == "Tr":
            self._require_arity(name, args, 1)
            value = args[0]
            if not _is_array(value):
                raise ValueError("Tr[...] requires a matrix.")
            return np.trace(value)

        if name == "Matrix":
            if not args:
                raise ValueError("Matrix[...] requires at least one factor.")

            result = args[0]
            if not _is_array(result):
                raise ValueError(
                    "Matrix[...] factors must be matrices."
                )

            for factor in args[1:]:
                if not _is_array(factor):
                    raise ValueError(
                        "Matrix[...] factors must be matrices."
                    )
                result = result @ factor

            return result

        # gen/heavy only occur in the ignored display-index bracket.  Reaching
        # them here would mean the export syntax changed in a meaningful way.
        raise ValueError(
            f"Unsupported RGBeta function head {name!r}."
        )

    @staticmethod
    def _require_arity(
        name: str,
        args: list[Value],
        count: int,
    ) -> None:
        if len(args) != count:
            raise ValueError(
                f"{name} expects {count} argument(s), got {len(args)}."
            )


def evaluate_inputform_expression(
    expression: str,
    environment: Mapping[str, Any],
) -> Value:
    """Evaluate one supported RGBeta ``InputForm`` expression."""

    return _Parser(str(expression), environment).parse()


def state_environment(
    state: T3UVState | SharedT3UVState,
) -> dict[str, Any]:
    """Return the RGBeta symbol -> numerical value environment."""

    state = state.validated()

    environment: dict[str, Any] = {
        "gY": state.sm.gY,
        "g2": state.sm.g2,
        "g3": state.sm.g3,
        "yu": state.sm.yu,
        "yd": state.sm.yd,
        "ye": state.sm.ye,
        "lambdaH": state.sm.lambdaH,
    }

    if isinstance(state, SharedT3UVState):
        environment.update(
            {
                "h": state.h,
                "MF": state.MF,
                "mSSq": state.mSSq,
                "lambdaS": state.lambdaS,
                "lambda3": state.lambda3,
                "lambda4": state.lambda4,
                "lambda5": state.lambda5,
            }
        )
        return environment

    environment.update(
        {
            "y1": state.y1,
            "y2": state.y2,
            "MF": state.MF,
            "mS1Sq": state.mS1Sq,
            "mS2Sq": state.mS2Sq,
            "lambdaS1": state.lambdaS1,
            "lambdaS2": state.lambdaS2,
            "lambdaH1": state.lambdaH1,
            "lambdaH2": state.lambdaH2,
            "lambda12": state.lambda12,
            "lambdaT3": state.lambdaT3,
        }
    )

    for name in (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross",
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    ):
        value = getattr(state, name)
        if value is not None:
            environment[name] = value

    return environment


def _validate_payload_metadata(
    payload: Mapping[str, Any],
    state: T3UVState | SharedT3UVState,
) -> None:
    metadata = payload.get("metadata", {})
    if not isinstance(metadata, Mapping):
        raise ValueError("RGBeta payload metadata must be an object.")

    rep = state.representation.validated()

    # Historical ordinary-T3 RGBeta exports did not include an explicit
    # SharedScalar field.  For those payloads, absence means the ordinary
    # split-scalar branch (False).  Shared-scalar exports must still declare
    # SharedScalar=True explicitly.
    actual_shared = bool(metadata.get("SharedScalar", False))

    expected = {
        "SharedScalar": rep.shared_scalar,
        "dS1": rep.d_s1,
        "dS2": rep.d_s2,
        "dF": rep.d_f,
        "alpha": rep.alpha,
    }

    actual_values = {
        "SharedScalar": actual_shared,
        "dS1": metadata.get("dS1"),
        "dS2": metadata.get("dS2"),
        "dF": metadata.get("dF"),
        "alpha": metadata.get("alpha"),
    }

    mismatches: list[str] = []
    for key, expected_value in expected.items():
        actual = actual_values[key]
        if actual != expected_value:
            mismatches.append(
                f"{key}: payload={actual!r}, state={expected_value!r}"
            )

    if mismatches:
        raise ValueError(
            "RGBeta payload does not match the numerical state: "
            + "; ".join(mismatches)
        )

    convention = metadata.get("ReportBetaConvention")
    expected_convention = (
        "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
    )
    if convention != expected_convention:
        raise ValueError(
            "Unsupported RGBeta report-beta convention. "
            f"Expected {expected_convention!r}, got {convention!r}."
        )


def evaluate_rgbeta_payload(
    payload: Mapping[str, Any],
    state: T3UVState | SharedT3UVState,
) -> dict[str, Value]:
    """Evaluate all UV ``report_betas`` for one numerical state."""

    if payload.get("status") != "Success":
        raise ValueError("RGBeta payload status is not Success.")

    state = state.validated()
    _validate_payload_metadata(payload, state)

    report_betas = payload.get("report_betas")
    if not isinstance(report_betas, Mapping):
        raise ValueError(
            "RGBeta payload does not contain a report_betas object."
        )

    environment = state_environment(state)

    evaluated = {
        str(name): evaluate_inputform_expression(
            str(expression),
            environment,
        )
        for name, expression in report_betas.items()
    }

    _, layout = pack_uv_state(state)
    validate_beta_keys(evaluated, layout)

    return evaluated


def load_rgbeta_payload(path: Path) -> dict[str, Any]:
    """Load one successful RGBeta JSON file."""

    path = Path(path)
    payload = json.loads(
        path.read_text(encoding="utf-8-sig")
    )

    if not isinstance(payload, dict):
        raise ValueError("RGBeta JSON top level must be an object.")

    return payload


def derivative_from_rgbeta_payload(
    payload: Mapping[str, Any],
    state: T3UVState | SharedT3UVState,
) -> tuple[np.ndarray, StateVectorLayout]:
    """Return ``dy/dln(mu)`` and its state-vector layout."""

    state = state.validated()
    _, layout = pack_uv_state(state)
    evaluated = evaluate_rgbeta_payload(payload, state)

    derivative = beta_values_to_derivative(
        evaluated,
        layout,
        divide_by_loop_factor=True,
    )

    return derivative, layout
