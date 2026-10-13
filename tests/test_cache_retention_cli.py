import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from zstarview import config
from zstarview.gui import viewer


@pytest.mark.parametrize("failure", [None, "read", "write", "mkdir"])
def test_cache_retention_command_reports_save_result(
    tmp_path, monkeypatch, capsys, failure
) -> None:
    cfg = tmp_path / "settings" / "config.json"
    monkeypatch.setattr(config, "_config_file", cfg)
    monkeypatch.setattr(
        viewer, "parse_args", lambda: SimpleNamespace(set_cache_retention_hours=48)
    )

    def unexpected_gui(*_args, **_kwargs):
        pytest.fail("The configuration command must not start the GUI")

    monkeypatch.setattr("zstarview.splash.setup_app", unexpected_gui)
    if failure == "read":
        cfg.mkdir(parents=True)
    elif failure == "mkdir":
        original_mkdir = Path.mkdir

        def denied_mkdir(path, *args, **kwargs):
            if path == cfg.parent:
                raise PermissionError("blocked")
            return original_mkdir(path, *args, **kwargs)

        monkeypatch.setattr(Path, "mkdir", denied_mkdir)
    elif failure == "write":
        original_write = Path.write_text

        def denied_write(path, *args, **kwargs):
            if path == cfg:
                raise PermissionError("blocked")
            return original_write(path, *args, **kwargs)

        monkeypatch.setattr(Path, "write_text", denied_write)

    with pytest.raises(SystemExit) as result:
        viewer.main()

    output = capsys.readouterr()
    if failure is None:
        assert result.value.code == 0
        assert output.out == "Cache history retention set to 48 hours.\n"
        assert output.err == ""
        assert json.loads(cfg.read_text(encoding="utf-8"))["cache_retention_hours"] == 48
    else:
        assert result.value.code == 1
        assert output.out == ""
        assert "Failed to save cache history retention settings." in output.err
