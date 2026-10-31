#!/usr/bin/env python3
"""Read-only diagnosis of indexed heavy-mass contractions in stage-1 WXF.

This does not run Mathematica, change model inputs, or implement a Wilson
coefficient. It parses a restricted subset of uncompressed WXF and fails closed
on unsupported tokens. It checks the top-level *tree* EFT only, where each
summand is a monomial (no nested symbolic sums).
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
import struct
import sys

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = (ROOT / "Reports/output/workstream_a/extended_c5_hierarchy/hierarchy/"
           "A_p0/seeded/debug/fresh_kernel_stage_2/input_transition_tree.wxf")


class WXFError(ValueError):
    pass


@dataclass(frozen=True)
class Symbol:
    name: str


@dataclass(frozen=True)
class Function:
    head: object
    args: tuple


def symbol_name(value: object) -> str:
    return value.name.split("`")[-1] if isinstance(value, Symbol) else ""


class WXFReader:
    """Only basic Wolfram expression records, never an evaluator."""
    def __init__(self, data: bytes):
        if len(data) > 10_000_000:
            raise WXFError("WXF input exceeds 10 MB safety bound")
        self.data = data
        self.pos = 0
        self.nodes = 0

    def take(self, count: int = 1) -> bytes:
        if count < 0 or self.pos + count > len(self.data):
            raise WXFError("Truncated WXF")
        b = self.data[self.pos:self.pos + count]
        self.pos += count
        return b

    def varint(self) -> int:
        result = 0
        for shift in range(0, 70, 7):
            n = self.take()[0]
            result |= (n & 127) << shift
            if n < 128:
                return result
        raise WXFError("Excessively long varint")

    def string(self) -> str:
        n = self.varint()
        if n > 10_000_000:
            raise WXFError("Oversized WXF string")
        return self.take(n).decode("utf-8")

    def expr(self, depth: int = 0) -> object:
        if depth > 200:
            raise WXFError("Excessive expression nesting")
        self.nodes += 1
        if self.nodes > 500_000:
            raise WXFError("Too many expression nodes")
        token = self.take()[0]
        if token == ord('f'):
            n = self.varint()
            if n > 100_000:
                raise WXFError("Oversized function")
            head = self.expr(depth + 1)
            return Function(head, tuple(self.expr(depth + 1) for _ in range(n)))
        if token == ord('s'):
            return Symbol(self.string())
        if token == ord('S'):
            return self.string()
        if token in (ord('C'), ord('j'), ord('i'), ord('L')):
            n = {ord('C'): 1, ord('j'): 2, ord('i'): 4, ord('L'): 8}[token]
            return int.from_bytes(self.take(n), "little", signed=True)
        if token == ord('r'):
            return struct.unpack('<d', self.take(8))[0]
        if token == ord('I'):
            return int(self.string())
        if token == ord('R'):
            return self.string()
        if token in (ord('A'), ord('-'), ord(':')):
            if token == ord('A'):
                n = self.varint()
                return Function(Symbol("Association"), tuple(self.expr(depth + 1) for _ in range(n)))
            return Function(Symbol("Rule" if token == ord('-') else "RuleDelayed"),
                            (self.expr(depth + 1), self.expr(depth + 1)))
        raise WXFError(f"Unsupported WXF token 0x{token:02x}; no unsafe fallback")

    def read(self) -> object:
        if self.take(2) != b'8:':
            raise WXFError("Expected uncompressed WXF v8 header")
        result = self.expr()
        if self.pos != len(self.data):
            raise WXFError("Unexpected trailing WXF bytes")
        return result


def walk(value: object):
    if isinstance(value, Function):
        yield value
        for a in value.args:
            yield from walk(a)


def index_key(node: Function):
    if symbol_name(node.head) == "Index" and len(node.args) == 2:
        idx, rep = node.args
        if not isinstance(idx, Symbol):
            raise WXFError("Unexpected index label")
        return (idx.name, flatten_name(rep))
    return None


def flatten_name(value: object) -> str:
    if isinstance(value, Symbol):
        return value.name
    if isinstance(value, Function):
        return flatten_name(value.head) + "[" + ",".join(map(flatten_name, value.args)) + "]"
    return str(value)


def coupling_names(term: object) -> list[str]:
    result = set()
    for node in walk(term):
        if symbol_name(node.head) == 'Coupling' and node.args:
            result.add(symbol_name(node.args[0]))
    return sorted(result)


def audit_tree(path: Path) -> dict:
    root = WXFReader(path.read_bytes()).read()
    if not isinstance(root, Function) or symbol_name(root.head) != 'Plus':
        raise WXFError("Expected top-level Plus of standalone tree EFT monomials")
    violations = []
    for position, term in enumerate(root.args, start=1):
        # The initial 43-term tree input has flat Times monomials. Counting
        # symbols inside expanded Plus would be mathematically meaningless.
        if not isinstance(term, Function) or symbol_name(term.head) != "Times":
            continue
        if any(symbol_name(node.head) == 'Plus' for node in walk(term)):
            continue
        counts = Counter(k for node in walk(term) if (k := index_key(node)) is not None)
        bad = [(label, rep, n) for (label, rep), n in counts.items() if n > 2]
        if bad:
            violations.append({
                "term_position_1_based": position,
                "overused_indices": [
                    {"index_symbol": label, "index_type": rep, "occurrences": n}
                    for label, rep, n in bad
                ],
                "couplings": coupling_names(term),
                "includes_indexed_heavy_mass_inverse": any(
                    symbol_name(node.head) == 'Power' and len(node.args) == 2 and
                    node.args[1] == -1 and isinstance(node.args[0], Function) and
                    symbol_name(node.args[0].head) == "Coupling" and
                    symbol_name(node.args[0].args[0]) == "MF"
                    for node in walk(term)
                ),
            })
    return {
        "status": "Stage2TreeIndexAuditComplete",
        "input_wxf": str(path.resolve()),
        "tree_monomials": len(root.args),
        "monomials_with_more_than_two_identical_typed_indices": len(violations),
        "violations": violations,
        "malformed_tree_is_valid_physical_prediction": False,
        "independent_mass_hierarchy_fixed": False,
        "suggested_exact_flavor_object": "W_pq = sum_r yA_pr yB_qr / MF_r",
    }


def assert_a0(result: dict) -> None:
    if result["tree_monomials"] != 43 or len(result["violations"]) != 4:
        raise WXFError("A_p0 stage-2 tree changed: not the 43-term/4-violation uploaded witness")
    if {v['term_position_1_based'] for v in result['violations']} != {33, 34, 42, 43}:
        raise WXFError("A_p0 violation positions changed")
    if not all(v['includes_indexed_heavy_mass_inverse'] and
               any(x['index_type'].endswith('`NFlavor') and x['occurrences'] == 3
                   for x in v['overused_indices'])
               for v in result['violations']):
        raise WXFError("A_p0 offending MF[r] index witness changed")
    if not any("y1Conj" in v['couplings'] for v in result['violations']):
        raise WXFError("Expected extra-Yukawa-dependent violation")


def self_test() -> None:
    tiny = b'8:f\x01s\x04Plusf\x02s\x05Timess\x01as\x01b'
    parsed = WXFReader(tiny).read()
    if not isinstance(parsed, Function) or len(parsed.args) != 1:
        raise AssertionError('Simple WXF parser failed')
    for corrupt in (b'9:f', tiny[:-1], tiny + b'junk', b'8:\xc1'):
        try:
            WXFReader(corrupt).read()
        except (WXFError, ValueError):
            continue
        raise AssertionError('Invalid WXF accepted')
    print('PASS: restricted WXF parser, malformed-data rejection')


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--wxf', type=Path, default=DEFAULT)
    p.add_argument('--output-json', type=Path)
    p.add_argument('--assert-a0', action='store_true')
    p.add_argument('--self-test', action='store_true')
    args = p.parse_args()
    try:
        self_test()
        if args.self_test:
            return 0
        result = audit_tree(args.wxf)
        if args.assert_a0:
            assert_a0(result)
        print(f"STAGE2 WXF: {result['tree_monomials']} tree monomials; "
              f"{len(result['violations'])} violate typed-index multiplicity")
        for v in result['violations']:
            print(f"  term {v['term_position_1_based']}: {', '.join(v['couplings'])}; "
                  f"indexed_MF_inverse={v['includes_indexed_heavy_mass_inverse']}")
        if args.output_json:
            args.output_json.parent.mkdir(parents=True, exist_ok=True)
            args.output_json.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
            print('Report:', args.output_json)
        return 0
    except (OSError, WXFError, ValueError) as exc:
        print('INDEX AUDIT BLOCKED: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
