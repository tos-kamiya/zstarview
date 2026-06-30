from __future__ import annotations

import logging

from zstarview import logging_utils


def test_setup_root_logger_can_skip_file_handler(monkeypatch) -> None:
    monkeypatch.setenv("ZSTARVIEW_DISABLE_FILE_LOGGING", "1")

    def fail_file_handler(*args, **kwargs):  # noqa: ANN001, ANN202
        raise AssertionError("FileHandler should not be created when disabled")

    monkeypatch.setattr(logging, "FileHandler", fail_file_handler)

    root_logger = logging_utils.setup_root_logger()

    assert root_logger is logging.getLogger()
