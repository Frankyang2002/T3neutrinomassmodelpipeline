#!/usr/bin/env python3
"""Encoding-safe entry point for the existing T3 full-general batch runner.

Uses the original runner unmodified. Only subprocess capture is adapted: Wolfram
can emit non-UTF-8 bytes (e.g. front-end box data) that crash CPython's
text-mode reader thread under -X utf8. Binary capture and backslashreplace
preserve diagnostics without fabricating a successful physical calculation.

Run from repository root:
    python -X utf8 -m scripts.RunT3ThesisFinalOneBatchSafe --all --run-id ... --debug
"""
from __future__ import annotations

import subprocess as _native_subprocess
from scripts import RunT3ThesisFinalOneBatch as _original


def _decode_output(value: bytes | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="backslashreplace")
    return value


class _BinaryCaptureAdapter:
    """Provide CompletedProcess text fields without text-mode reader threads."""
    TimeoutExpired = _native_subprocess.TimeoutExpired
    CalledProcessError = _native_subprocess.CalledProcessError
    DEVNULL = _native_subprocess.DEVNULL
    PIPE = _native_subprocess.PIPE

    @staticmethod
    def run(*args, **kwargs):
        original_text = kwargs.pop("text", False)
        universal = kwargs.pop("universal_newlines", False)
        enc = kwargs.pop("encoding", None)
        kwargs.pop("errors", None)
        expects_text = bool(original_text or universal or enc)
        # Keep native stdout and stderr as bytes until the entire child exits.
        result = _native_subprocess.run(*args, **kwargs)
        if not expects_text:
            return result
        return _native_subprocess.CompletedProcess(
            args=result.args,
            returncode=result.returncode,
            stdout=_decode_output(result.stdout),
            stderr=_decode_output(result.stderr),
        )


# Only the imported runner's local module reference is replaced. The Python
# standard-library subprocess module is not monkey-patched globally.
_original.subprocess = _BinaryCaptureAdapter


def main(argv: list[str] | None = None) -> int:
    return _original.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
