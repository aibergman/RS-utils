"""Explicit, reusable filesystem and workspace primitives.

These helpers never change the caller's process working directory.  Legacy
fixed-filename programs can be run by callers with ``subprocess`` and
``cwd=workspace.path`` while the rest of the application keeps explicit paths.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Union


PathLike = Union[str, os.PathLike[str]]


def as_path(path: PathLike) -> Path:
    """Return a ``Path`` without consulting or changing the process CWD."""

    return Path(path)


class Workspace:
    """An owned temporary directory that is cleaned up on context exit."""

    def __init__(self, path: Path, *, keep: bool = False, owned: bool = True):
        self.path = Path(path)
        self.keep = keep
        self._owned = owned
        self._cleaned = False

    @classmethod
    def create(
        cls,
        parent: PathLike | None = None,
        *,
        prefix: str = "rsutil-",
        same_filesystem_as: PathLike | None = None,
        keep: bool = False,
    ) -> "Workspace":
        """Create an isolated directory, optionally beside a destination.

        ``same_filesystem_as`` is useful for staging a candidate before an
        atomic replacement.  The temporary directory is created in the
        destination's parent directory, so ``os.replace`` can remain atomic.
        """

        if parent is not None and same_filesystem_as is not None:
            raise ValueError("pass either parent or same_filesystem_as, not both")
        if same_filesystem_as is not None:
            anchor = Path(same_filesystem_as)
            parent_path = anchor.parent
        elif parent is not None:
            parent_path = Path(parent)
        else:
            parent_path = None
        if parent_path is not None:
            parent_path.mkdir(parents=True, exist_ok=True)
        path = Path(tempfile.mkdtemp(prefix=prefix, dir=str(parent_path) if parent_path else None))
        return cls(path, keep=keep)

    def __enter__(self) -> Path:
        self.path.mkdir(parents=True, exist_ok=True)
        return self.path

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.cleanup()

    def cleanup(self) -> None:
        """Remove this workspace, including after an exception.

        Only the directory created by ``Workspace.create`` is removable by
        this object.  Cleanup is idempotent and can be disabled with
        ``keep=True`` for debugging.
        """

        if self._cleaned or self.keep or not self._owned:
            return
        shutil.rmtree(self.path)
        self._cleaned = True

    def child(self, name: str) -> Path:
        """Return a resolved child path and reject lexical/symlink escapes."""

        candidate = self.path / name
        root = self.path.resolve()
        resolved = candidate.resolve(strict=False)
        try:
            resolved.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"workspace child escapes workspace: {name!r}") from exc
        return resolved


@contextmanager
def isolated_workspace(
    parent: PathLike | None = None,
    *,
    prefix: str = "rsutil-",
    same_filesystem_as: PathLike | None = None,
) -> Iterator[Path]:
    """Yield an isolated temporary workspace and always clean it up."""

    workspace = Workspace.create(
        parent, prefix=prefix, same_filesystem_as=same_filesystem_as
    )
    try:
        yield workspace.path
    finally:
        workspace.cleanup()


def staging_workspace(destination: PathLike, *, prefix: str = "rsutil-stage-"):
    """Return a context manager for a same-filesystem staging directory."""

    return isolated_workspace(prefix=prefix, same_filesystem_as=destination)


def stage_file(source: PathLike, workspace: PathLike, *, name: str | None = None) -> Path:
    """Copy one explicit source file into a workspace and return its path."""

    source_path = Path(source)
    workspace_path = Path(workspace)
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    workspace_path.mkdir(parents=True, exist_ok=True)
    target = workspace_path / (name if name is not None else source_path.name)
    if target.parent != workspace_path:
        raise ValueError("staged file name must not contain a directory")
    shutil.copy2(source_path, target)
    return target


def stage_files(sources: Iterable[PathLike], workspace: PathLike) -> tuple[Path, ...]:
    """Copy explicit source files into a workspace, preserving file names."""

    return tuple(stage_file(source, workspace) for source in sources)


def atomic_replace(source: PathLike, destination: PathLike, *, overwrite: bool = True) -> Path:
    """Atomically replace a destination with a staged file on the same filesystem."""

    source_path = Path(source)
    destination_path = Path(destination)
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    if not overwrite and destination_path.exists():
        raise FileExistsError(destination_path)
    os.replace(source_path, destination_path)
    return destination_path


temporary_workspace = isolated_workspace
