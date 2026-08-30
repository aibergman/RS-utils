# Test dependency and fixtures

Install the documented development extra and run the active regression suite
with:

```bash
python -m pip install -e '.[dev]'
python -m pytest
```

The tests copy only the selected files from `examples/` into pytest temporary
directories. They never modify the repository examples in place.
