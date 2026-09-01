# Test dependency and fixtures

Install the documented test extra and run the active regression suite
with:

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```

The tests copy only the selected files from `examples/` into pytest temporary
directories. They never modify the repository examples in place.
