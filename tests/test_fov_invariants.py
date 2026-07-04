from zstarview.types import ViewerData, ViewProjection


def test_viewer_data_clamps_content_fov_to_edge_fov() -> None:
    viewer = ViewerData(
        location=(35.0, 139.0),
        timezone_name="UTC",
        city_name="Tokyo",
        edge_fov_deg=110.0,
        content_fov_deg=95.0,
    )

    assert viewer.edge_fov_deg == 110.0
    assert viewer.content_fov_deg == 110.0


def test_view_projection_clamps_content_fov_to_edge_fov() -> None:
    projection = ViewProjection(
        view_center=(45.0, 180.0),
        edge_fov_deg=120.0,
        content_fov_deg=100.0,
    )

    assert projection.edge_fov_deg == 120.0
    assert projection.content_fov_deg == 120.0
