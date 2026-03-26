from .location_resolver import viewpoints as _impl

Viewpoint = _impl.Viewpoint
ascii_fallback_name = _impl.ascii_fallback_name
build_viewpoint = _impl.build_viewpoint
find_exact_viewpoint_matches = _impl.find_exact_viewpoint_matches
list_viewpoint_all_names = _impl.list_viewpoint_all_names
list_viewpoint_primary_names = _impl.list_viewpoint_primary_names
load_viewpoints = _impl.load_viewpoints
looks_like_qid_placeholder = _impl.looks_like_qid_placeholder
normalize_viewpoint_name = _impl.normalize_viewpoint_name
prefixed_viewpoint_name = _impl.prefixed_viewpoint_name
resolve_viewpoint = _impl.resolve_viewpoint
split_prefixed_viewpoint = _impl.split_prefixed_viewpoint
viewpoint_to_dict = _impl.viewpoint_to_dict

__all__ = [
    "Viewpoint",
    "ascii_fallback_name",
    "build_viewpoint",
    "find_exact_viewpoint_matches",
    "list_viewpoint_all_names",
    "list_viewpoint_primary_names",
    "load_viewpoints",
    "looks_like_qid_placeholder",
    "normalize_viewpoint_name",
    "prefixed_viewpoint_name",
    "resolve_viewpoint",
    "split_prefixed_viewpoint",
    "viewpoint_to_dict",
]
