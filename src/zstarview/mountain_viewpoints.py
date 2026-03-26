from .location_resolver import mountains as _impl

MountainViewpoint = _impl.MountainViewpoint
list_mountain_all_names = _impl.list_mountain_all_names
list_mountain_primary_names = _impl.list_mountain_primary_names
load_mountain_viewpoints = _impl.load_mountain_viewpoints
mountain_viewpoint_to_dict = _impl.mountain_viewpoint_to_dict
resolve_mountain_viewpoint = _impl.resolve_mountain_viewpoint

__all__ = [
    "MountainViewpoint",
    "list_mountain_all_names",
    "list_mountain_primary_names",
    "load_mountain_viewpoints",
    "mountain_viewpoint_to_dict",
    "resolve_mountain_viewpoint",
]
