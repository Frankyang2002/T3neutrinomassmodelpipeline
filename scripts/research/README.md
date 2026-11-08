# Research-only scripts

This directory is separate from the supported pipeline entry points:
`pipeline.py` and `full16_symbolic.py` in the project root. The command index
at [`scripts/INDEX.md`](../INDEX.md) classifies the remaining scripts by
scientific responsibility without changing their established import paths.

## Experimental programs

The optional `experimental/` package is for independent, standalone
investigations. A file is moved here only after a local audit finds no
retained references to its former pathname, no use of `__file__`, and no
obvious dynamic/relative import dependency. These checks do **not** prove
complete runtime equivalence, so experimental relocation is opt-in.

Run a relocated program from the project root using Python's module interface,
for example:

```powershell
python -X utf8 -m scripts.research.experimental.ExampleProgram --help
```

`ExampleProgram` is illustrative. Use an actual filename reported in the
reviewed move manifest. The migration preserves original source bytes; the
`.experimental.py` suffix becomes a normal `.py` filename because the
experimental status is expressed by its directory.
