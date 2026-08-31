"""Release packaging and version contract checks."""

import rsutil


def test_package_version_is_release_candidate() -> None:
    assert rsutil.__version__ == "0.11.0rc1"
