import pytest

from zstarview.cli.args import parse_args


def test_content_fov_default_is_current_overscan() -> None:
    args = parse_args([])
    assert args.content_fov_deg == 115.0


def test_content_fov_accepts_valid_overscan_value() -> None:
    args = parse_args(["--content-fov-deg", "110"])
    assert args.content_fov_deg == 110.0


@pytest.mark.parametrize("value", ["89.9", "127.1"])
def test_content_fov_rejects_out_of_range_values(value: str) -> None:
    with pytest.raises(SystemExit):
        parse_args(["--content-fov-deg", value])
