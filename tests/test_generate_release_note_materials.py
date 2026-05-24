from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "generate_release_note_materials.py"
SPEC = importlib.util.spec_from_file_location("generate_release_note_materials", SCRIPT_PATH)
assert SPEC is not None
assert SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_parse_version_from_text() -> None:
    text = '# SPDX-License-Identifier: MIT\n__version__ = "1.2.3"\n'
    assert MODULE._parse_version_from_text(text) == "1.2.3"


def test_is_bump_only_subject() -> None:
    assert MODULE._is_bump_only_subject("chore: bump version to 1.2.3")
    assert MODULE._is_bump_only_subject("fix: bump version to 1.2.3")
    assert MODULE._is_bump_only_subject("build: bump version to 1.2.3")
    assert MODULE._is_bump_only_subject("fix: revert version bump")
    assert not MODULE._is_bump_only_subject("feat: add release-note generator")


def test_format_commit_list_uses_short_hash() -> None:
    commit = MODULE.CommitInfo(
        commit_hash="1234567890abcdef",
        committed_at="2026-05-25T06:00:00+09:00",
        subject="feat: add release-note generator",
    )
    assert MODULE._format_commit_list([commit]) == [
        "- 12345678 2026-05-25T06:00:00+09:00 feat: add release-note generator"
    ]
