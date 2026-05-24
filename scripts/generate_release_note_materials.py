#!/usr/bin/env python3
"""Generate release-note source materials from git history.

This script treats the first commit that changes `__version__` to a released
version as the release commit. It then gathers the commit subjects since the
previous released version and the `pyproject.toml` diff across that release
window.

Example:
  python scripts/generate_release_note_materials.py
  python scripts/generate_release_note_materials.py --output /tmp/release-notes.md
  python scripts/generate_release_note_materials.py --include-bump-commits
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_VERSION_FILE = Path("src/zstarview/__about__.py")
DEFAULT_VERSIONS_FILE = Path("releases.txt")
DEFAULT_PYPROJECT_FILE = Path("pyproject.toml")

VERSION_RE = re.compile(r'__version__\s*=\s*"([^"]+)"')
BUMP_SUBJECT_RE = re.compile(
    r"^(?:chore|fix|build): bump version to \S+$|^fix: revert version bump$"
)


@dataclass(frozen=True)
class CommitInfo:
    commit_hash: str
    committed_at: str
    subject: str

    @property
    def short_hash(self) -> str:
        return self.commit_hash[:8]


@dataclass(frozen=True)
class ReleaseInfo:
    version: str
    release_commit: CommitInfo
    release_body: str


def _run_git(args: Sequence[str]) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"git {' '.join(args)} failed")
    return proc.stdout


def _parse_version_from_text(text: str) -> str:
    match = VERSION_RE.search(text)
    if not match:
        raise ValueError("Could not find __version__ assignment")
    return match.group(1)


def _parse_commit_line(line: str) -> CommitInfo:
    commit_hash, committed_at, subject = line.split("\t", 2)
    return CommitInfo(
        commit_hash=commit_hash,
        committed_at=committed_at,
        subject=subject,
    )


def _read_released_versions(path: Path) -> list[str]:
    versions = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    return [version for version in versions if version]


def _get_commit_info(commit_hash: str) -> CommitInfo:
    line = _run_git(["show", "-s", "--format=%H%x09%cI%x09%s", commit_hash]).strip()
    return _parse_commit_line(line)


def _get_commit_body(commit_hash: str) -> str:
    body = _run_git(["show", "-s", "--format=%b", commit_hash])
    return body.strip()


def _get_file_text_at_commit(commit_hash: str, path: Path) -> str:
    return _run_git(["show", f"{commit_hash}:{path.as_posix()}"])


def _collect_version_transitions(version_file: Path) -> list[ReleaseInfo]:
    log_output = _run_git(["log", "--follow", "--reverse", "--format=%H", "--", version_file.as_posix()])
    commit_hashes = [line.strip() for line in log_output.splitlines() if line.strip()]
    transitions: list[ReleaseInfo] = []
    current_version: str | None = None
    for commit_hash in commit_hashes:
        text = _get_file_text_at_commit(commit_hash, version_file)
        version = _parse_version_from_text(text)
        if version == current_version:
            continue
        commit_info = _get_commit_info(commit_hash)
        transitions.append(
            ReleaseInfo(
                version=version,
                release_commit=commit_info,
                release_body=_get_commit_body(commit_hash),
            )
        )
        current_version = version
    return transitions


def _is_bump_only_subject(subject: str) -> bool:
    return bool(BUMP_SUBJECT_RE.match(subject.strip()))


def _get_commits_in_range(older_hash: str | None, newer_hash: str) -> list[CommitInfo]:
    if older_hash is None:
        args = ["log", "--reverse", "--format=%H%x09%cI%x09%s", newer_hash]
    else:
        args = ["log", "--reverse", "--format=%H%x09%cI%x09%s", f"{older_hash}..{newer_hash}"]
    output = _run_git(args)
    return [_parse_commit_line(line) for line in output.splitlines() if line.strip()]


def _get_pyproject_diff(older_hash: str | None, newer_hash: str, pyproject_file: Path) -> str:
    if older_hash is None:
        return _get_file_text_at_commit(newer_hash, pyproject_file).strip()
    diff_text = _run_git(["diff", older_hash, newer_hash, "--", pyproject_file.as_posix()]).strip()
    return diff_text


def _format_commit_list(commits: Sequence[CommitInfo]) -> list[str]:
    return [f"- {commit.short_hash} {commit.committed_at} {commit.subject}" for commit in commits]


def _format_release_section(
    release: ReleaseInfo,
    previous_release: ReleaseInfo | None,
    commits_in_range: Sequence[CommitInfo],
    pyproject_diff: str,
    *,
    include_bump_commits: bool,
) -> str:
    notable_commits = list(commits_in_range)
    if not include_bump_commits:
        notable_commits = [commit for commit in commits_in_range if not _is_bump_only_subject(commit.subject)]
    bump_commits = [commit for commit in commits_in_range if _is_bump_only_subject(commit.subject)]

    lines = [f"## {release.version}", ""]
    lines.append(f"Release commit: `{release.release_commit.commit_hash}`")
    lines.append(f"Release date: `{release.release_commit.committed_at}`")
    lines.append(f"Release subject: {release.release_commit.subject}")
    if previous_release is None:
        lines.append("Previous released version: `(none in input list)`")
    else:
        lines.append(
            "Previous released version: "
            f"`{previous_release.version}` (`{previous_release.release_commit.commit_hash}`)"
        )
    lines.append(
        "Commits in release window: "
        f"{len(commits_in_range)} total, {len(notable_commits)} shown"
    )
    if not include_bump_commits:
        lines.append(f"Skipped bump-only commits: {len(bump_commits)}")
    lines.append("")

    if release.release_body:
        lines.append("Release commit body:")
        lines.append("```text")
        lines.append(release.release_body)
        lines.append("```")
        lines.append("")

    lines.append("Commit subjects in release window:")
    if notable_commits:
        lines.extend(_format_commit_list(notable_commits))
    else:
        lines.append("- (none)")
    lines.append("")

    if not include_bump_commits and bump_commits:
        lines.append("Version-bump commits in release window:")
        lines.extend(_format_commit_list(bump_commits))
        lines.append("")

    if previous_release is None:
        lines.append("Initial `pyproject.toml` snapshot at release commit:")
        fence = "toml"
    else:
        lines.append("`pyproject.toml` diff from previous released version:")
        fence = "diff"
    lines.append(f"```{fence}")
    lines.append(pyproject_diff if pyproject_diff else "(no pyproject.toml changes)")
    lines.append("```")
    lines.append("")
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate release-note source materials from version bumps and git history."
    )
    parser.add_argument(
        "--versions-file",
        type=Path,
        default=DEFAULT_VERSIONS_FILE,
        help=f"Released versions list file (default: {DEFAULT_VERSIONS_FILE}).",
    )
    parser.add_argument(
        "--version-file",
        type=Path,
        default=DEFAULT_VERSION_FILE,
        help=f"Path to the version source file (default: {DEFAULT_VERSION_FILE}).",
    )
    parser.add_argument(
        "--pyproject-file",
        type=Path,
        default=DEFAULT_PYPROJECT_FILE,
        help=f"Path to pyproject.toml (default: {DEFAULT_PYPROJECT_FILE}).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output path. Defaults to stdout.",
    )
    parser.add_argument(
        "--include-bump-commits",
        action="store_true",
        help="Include version-bump-only commits in the main commit list.",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()

    released_versions = _read_released_versions(args.versions_file)
    transitions = _collect_version_transitions(args.version_file)
    transition_by_version = {release.version: release for release in transitions}

    missing_versions = [version for version in released_versions if version not in transition_by_version]
    if missing_versions:
        missing = ", ".join(missing_versions)
        raise SystemExit(f"Versions not found in git history for {args.version_file}: {missing}")

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    sections = [
        "# zstarview release-note materials",
        "",
        f"Generated at: `{generated_at}`",
        f"Versions file: `{args.versions_file.as_posix()}`",
        f"Version source: `{args.version_file.as_posix()}`",
        f"Pyproject path: `{args.pyproject_file.as_posix()}`",
        "",
        (
            "Each section treats the first commit that changed `__version__` "
            "to the released version as the release commit."
        ),
        "",
    ]

    selected_releases = [transition_by_version[version] for version in released_versions]
    for index, release in enumerate(selected_releases):
        previous_release = selected_releases[index + 1] if index + 1 < len(selected_releases) else None
        older_hash = previous_release.release_commit.commit_hash if previous_release else None
        newer_hash = release.release_commit.commit_hash
        commits_in_range = _get_commits_in_range(older_hash, newer_hash)
        pyproject_diff = _get_pyproject_diff(older_hash, newer_hash, args.pyproject_file)
        sections.append(
            _format_release_section(
                release,
                previous_release,
                commits_in_range,
                pyproject_diff,
                include_bump_commits=args.include_bump_commits,
            )
        )

    output_text = "\n".join(sections).rstrip() + "\n"
    if args.output is None:
        sys.stdout.write(output_text)
    else:
        args.output.write_text(output_text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
