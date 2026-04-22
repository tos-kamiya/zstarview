import pytest

from zstarview.cli.args import parse_args


def test_content_fov_default_is_current_overscan() -> None:
    args = parse_args([])
    assert args.content_fov_deg == 110.0
    assert args.edge_fov_deg == 95.0


def test_content_fov_accepts_valid_overscan_value() -> None:
    args = parse_args(["--content-fov-deg", "110"])
    assert args.content_fov_deg == 110.0


@pytest.mark.parametrize("value", ["89.9", "135.1"])
def test_content_fov_rejects_out_of_range_values(value: str) -> None:
    with pytest.raises(SystemExit):
        parse_args(["--content-fov-deg", value])


def test_edge_fov_rejects_values_larger_than_content_fov() -> None:
    with pytest.raises(SystemExit):
        parse_args(["--edge-fov-deg", "120", "--content-fov-deg", "110"])


def test_edge_fov_allows_values_up_to_content_fov() -> None:
    args = parse_args(["--edge-fov-deg", "110", "--content-fov-deg", "110"])
    assert args.edge_fov_deg == 110.0
