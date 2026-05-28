from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Literal

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QLineEdit,
    QFrame,
    QPushButton,
    QScrollArea,
    QRadioButton,
    QSpinBox,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from ..cli.args import (
    _parse_cloud_stripe,
    _parse_window_geometry,
)
from .place_search_dialog import PlaceSearchDialog
from .famous_star_shortcuts import build_place_search_jump_targets
from ..location_resolver import ResolvedLocation, resolve_launch_location
from ..location_resolver.place_search import search_place_candidates
from ..paths import (
    OVERLAY_FONT_SIZE_MAX,
    OVERLAY_FONT_SIZE_MIN,
    THEME_PRESET_NAMES,
)
from .launch_profile import default_gui_launch_profile
from .worker_pool import submit_gui_work

TriBool = Literal["default", "true", "false"]
_FLOAT_DEFAULTS: dict[str, float] = {
    "observer_height_m": 1.7,
}
_INT_DEFAULTS: dict[str, int] = {}


@dataclass(frozen=True, slots=True)
class _FieldSpec:
    key: str
    label: str
    kind: Literal["text", "float", "int", "bool", "choice", "note"]
    tab: str
    choices: tuple[str, ...] = ()
    decimals: int = 3
    minimum: float = 0.0
    maximum: float = 1.0
    step: float = 1.0
    placeholder: str = ""
    note: str = ""


class _CollapsibleSection(QWidget):
    def __init__(self, title: str, *, expanded: bool = True, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._button = QToolButton(self)
        self._button.setText(title)
        self._button.setCheckable(True)
        self._button.setChecked(expanded)
        self._button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._button.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self._button.toggled.connect(self._set_expanded)

        self._content = QWidget(self)
        self._content_layout = QFormLayout(self._content)
        self._content_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
        self._content_layout.setFormAlignment(
            Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
        )
        self._content_layout.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        layout.addWidget(self._button)
        layout.addWidget(self._content)
        self._set_expanded(expanded)

    @property
    def form_layout(self) -> QFormLayout:
        return self._content_layout

    def is_expanded(self) -> bool:
        return bool(self._button.isChecked())

    def _set_expanded(self, expanded: bool) -> None:
        self._button.setArrowType(
            Qt.ArrowType.DownArrow if expanded else Qt.ArrowType.RightArrow
        )
        self._content.setVisible(expanded)


def _as_text(value: Any, *, key: str) -> str:
    if value is None:
        if key == "window_geometry":
            return "restore"
        return ""
    if key == "cloud_stripe" and isinstance(value, (tuple, list)) and len(value) == 3:
        mode, count, width = value
        return f"{mode},{count},{width}"
    if key == "window_geometry" and isinstance(value, (tuple, list)) and len(value) == 4:
        x, y, width, height = value
        return f"{x},{y},{width},{height}"
    return str(value)


def _tri_bool_to_value(text: str) -> bool | None:
    lowered = text.strip().lower()
    if lowered == "default":
        return None
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    raise ValueError(f"Invalid tri-state boolean: {text!r}")


def _coerce_float_value(key: str, value: Any, fallback: float) -> float:
    if value is None:
        return float(_FLOAT_DEFAULTS.get(key, fallback))
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(_FLOAT_DEFAULTS.get(key, fallback))


def _coerce_int_value(key: str, value: Any, fallback: int) -> int:
    if value is None:
        return int(_INT_DEFAULTS.get(key, fallback))
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(_INT_DEFAULTS.get(key, fallback))


class StartupDialog(QDialog):
    city_auto_finished = Signal(int, object, str)

    def __init__(
        self,
        profile: dict[str, Any] | None = None,
        *,
        auto_location_resolver: Callable[[], ResolvedLocation] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("zstarview-gui startup")
        self.setModal(True)
        self.resize(480, 456)

        self._defaults = default_gui_launch_profile()
        self._base_profile = dict(self._defaults)
        if profile:
            self._base_profile.update(profile)

        self._widgets: dict[str, Any] = {}
        self._tab_layouts: dict[str, QFormLayout] = {}
        self._location_group_widgets: dict[str, list[QWidget]] = {"city": [], "place": []}
        self._time_group_layouts: dict[str, QFormLayout] = {}
        self._time_group_widgets: dict[str, list[QWidget]] = {}
        self._overlay_layouts: dict[str, QFormLayout] = {}
        self._overlay_sections: dict[str, _CollapsibleSection] = {}
        self._overlay_section_by_key: dict[str, str] = {}
        self._city_auto_request_id = 0
        self._auto_location_resolver = auto_location_resolver or _resolve_auto_location_for_dialog
        self.city_auto_finished.connect(self._on_city_auto_finished)
        self._city_auto_button: QPushButton | None = None
        self._city_search_button: QPushButton | None = None
        self._place_selected_payload: dict[str, Any] | None = None
        self._location_city_radio: QRadioButton | None = None
        self._location_place_radio: QRadioButton | None = None
        self._location_mode_button_group = QButtonGroup(self)
        self._time_source_button_group = QButtonGroup(self)
        self._time_current_radio: QRadioButton | None = None
        self._time_relative_radio: QRadioButton | None = None
        self._time_absolute_radio: QRadioButton | None = None
        self._view_center_az_widget: QDoubleSpinBox | None = None
        self._view_center_alt_hint_label: QLabel | None = None
        self._view_center_az_buttons: dict[str, QPushButton] = {}

        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(12, 12, 12, 12)

        self._tabs = QTabWidget(self)
        outer_layout.addWidget(self._tabs, 1)

        tab_order = (
            "Location",
            "View",
            "Time",
            "Stars",
            "Overlays",
            "General",
        )
        for tab_name in tab_order:
            tab_widget = QWidget(self)
            scroll_area = QScrollArea(self)
            scroll_area.setWidgetResizable(True)
            scroll_area.setFrameShape(QFrame.Shape.NoFrame)
            if tab_name == "Location":
                tab_layout = QFormLayout(tab_widget)
                tab_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
                tab_layout.setFormAlignment(
                    Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
                )
                tab_layout.setFieldGrowthPolicy(
                    QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
                )
                self._tab_layouts[tab_name] = tab_layout
                self._build_location_mode_selector(tab_widget, tab_layout)
            elif tab_name == "View":
                tab_layout = QFormLayout(tab_widget)
                tab_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
                tab_layout.setFormAlignment(
                    Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
                )
                tab_layout.setFieldGrowthPolicy(
                    QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
                )
                self._tab_layouts[tab_name] = tab_layout
            elif tab_name == "Time":
                tab_layout = QVBoxLayout(tab_widget)
                tab_layout.setContentsMargins(0, 0, 0, 0)
                tab_layout.setSpacing(12)
                self._build_time_tab(tab_widget, tab_layout)
            elif tab_name == "Overlays":
                tab_layout = QVBoxLayout(tab_widget)
                tab_layout.setContentsMargins(0, 0, 0, 0)
                tab_layout.setSpacing(10)
                self._build_overlay_tab(tab_widget, tab_layout)
            else:
                tab_layout = QFormLayout(tab_widget)
                tab_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
                tab_layout.setFormAlignment(
                    Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
                )
                tab_layout.setFieldGrowthPolicy(
                    QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
                )
                self._tab_layouts[tab_name] = tab_layout
            scroll_area.setWidget(tab_widget)
            self._tabs.addTab(scroll_area, tab_name)

        self._add_specs()
        self._restore_from_profile(self._base_profile)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self._reset_button = QPushButton("Reset to Default Values", self)
        self._reset_button.clicked.connect(self.reset_to_defaults)
        button_row.addWidget(self._reset_button)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        button_row.addWidget(buttons)
        outer_layout.addLayout(button_row)

    def _add_specs(self) -> None:
        specs = (
            _FieldSpec("city", "City", "text", "Location"),
            _FieldSpec("observer_height_m", "Observer height", "float", "Location", minimum=0.0, maximum=10000.0, step=0.1),
            _FieldSpec("use_building_top", "Use building top", "bool", "Location"),
            _FieldSpec("view_center_alt", "View center alt", "float", "View", minimum=-90.0, maximum=90.0, step=1.0),
            _FieldSpec("view_center_az", "View center az", "float", "View", minimum=0.0, maximum=360.0, step=1.0),
            _FieldSpec("edge_fov_deg", "Edge FOV", "float", "View", minimum=0.1, maximum=135.0, step=0.5),
            _FieldSpec("content_fov_deg", "Content FOV", "float", "View", minimum=90.0, maximum=135.0, step=0.5),
            _FieldSpec("hours", "Hours", "float", "Time", minimum=-9999.0, maximum=9999.0, step=0.5),
            _FieldSpec("days", "Days", "float", "Time", minimum=-9999.0, maximum=9999.0, step=0.5),
            _FieldSpec("datetime", "Date/time", "text", "Time"),
            _FieldSpec("timezone", "Timezone", "text", "Time"),
            _FieldSpec("sky_opacity", "Sky opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("sky_disc_style", "Sky disc style", "choice", "Overlays", choices=("smooth",)),
            _FieldSpec("sky_disc_altaz_rings", "Sky disc rings", "choice", "Overlays", choices=("off", "dimalt", "altaz")),
            _FieldSpec("sky_disc_altaz_rings_hover", "Hover rings", "choice", "Overlays", choices=("off", "dimalt", "altaz")),
            _FieldSpec("vmag_limit", "Vmag limit", "float", "Stars", minimum=0.0, maximum=20.0, step=0.1),
            _FieldSpec("vmag_brightness_multiplier", "Brightness multiplier", "float", "Stars", minimum=1.0, maximum=3.0, step=0.01),
            _FieldSpec("enlarge_moon", "Enlarge moon", "bool", "Stars"),
            _FieldSpec("bright_bodies", "Bright bodies", "choice", "Stars", choices=("outline", "fill")),
            _FieldSpec("star_base_radius", "Star base radius", "float", "Stars", minimum=0.0, maximum=20.0, step=0.1),
            _FieldSpec("expected_render_width", "Expected render width", "int", "Stars", minimum=1.0, maximum=10000.0, step=1.0),
            _FieldSpec("show_dso_initial", "DSO visibility", "choice", "Stars", choices=("default", "true", "false")),
            _FieldSpec("show_asterisms_initial", "Asterisms visibility", "choice", "Stars", choices=("default", "true", "false")),
            _FieldSpec("show_guidelines_initial", "Guidelines visibility", "choice", "Stars", choices=("default", "true", "false")),
            _FieldSpec("theme", "Theme", "choice", "General", choices=THEME_PRESET_NAMES),
            _FieldSpec("window_geometry", "Window geometry", "text", "General"),
            _FieldSpec("window_frame", "Window frame", "choice", "General", choices=("frameless", "window")),
            _FieldSpec("observation_info", "Observation info", "choice", "General", choices=("auto", "top", "bottom", "off")),
            _FieldSpec("visibility_boost", "Visibility boost", "float", "General", minimum=1.0, maximum=10.0, step=0.1),
            _FieldSpec("overlay_font_size", "Overlay font size", "float", "General", minimum=float(OVERLAY_FONT_SIZE_MIN), maximum=float(OVERLAY_FONT_SIZE_MAX), step=0.5),
            _FieldSpec("cloud_opacity", "Cloud opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("cloud_stripe", "Cloud stripe", "text", "Overlays"),
            _FieldSpec("cloud_missing_tint_opacity", "Cloud missing tint", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("aircraft_opacity", "Aircraft opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("satellite_opacity", "Satellite opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("terrain_horizon_opacity", "Terrain horizon opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.001),
            _FieldSpec("earth_guide_opacity", "Earth guide opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.001),
            _FieldSpec("night_light_opacity", "Night light opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("ground_tint_opacity", "Ground tint opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("urban_outline_opacity", "Urban outline opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("water_surface_opacity", "Water surface opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("urban_outline_feature_type", "Urban outline mode", "choice", "Overlays", choices=("both", "building")),
            _FieldSpec("urban_outline_radius_km", "Urban radius km", "float", "Overlays", minimum=0.0, maximum=1000.0, step=0.1),
            _FieldSpec("urban_outline_skyscraper_radius_km", "Skyscraper radius km", "float", "Overlays", minimum=0.0, maximum=1000.0, step=0.1),
            _FieldSpec("urban_outline_min_height_m", "Min building height", "float", "Overlays", minimum=0.0, maximum=100000.0, step=0.1),
            _FieldSpec("urban_outline_skyscraper_only", "Skyscraper only", "bool", "Overlays"),
        )
        for spec in specs:
            self._add_spec(spec)

    def _build_overlay_tab(self, tab_widget: QWidget, layout: QVBoxLayout) -> None:
        section_defs = (
            ("Sky", ("sky_opacity", "sky_disc_style", "sky_disc_altaz_rings", "sky_disc_altaz_rings_hover")),
            ("Clouds", ("cloud_opacity", "cloud_stripe", "cloud_missing_tint_opacity")),
            ("Aircraft and Satellites", ("aircraft_opacity", "satellite_opacity")),
            (
                "Ground and Guides",
                (
                    "terrain_horizon_opacity",
                    "earth_guide_opacity",
                    "night_light_opacity",
                    "ground_tint_opacity",
                    "water_surface_opacity",
                ),
            ),
            (
                "Urban Outline",
                (
                    "urban_outline_opacity",
                    "urban_outline_feature_type",
                    "urban_outline_radius_km",
                    "urban_outline_skyscraper_radius_km",
                    "urban_outline_min_height_m",
                    "urban_outline_skyscraper_only",
                ),
            ),
        )
        for title, keys in section_defs:
            section = _CollapsibleSection(title, parent=tab_widget)
            layout.addWidget(section)
            self._overlay_sections[title] = section
            self._overlay_layouts[title] = section.form_layout
            for key in keys:
                self._overlay_section_by_key[key] = title
        layout.addStretch(1)

    def _build_location_mode_selector(self, tab_widget: QWidget, layout: QFormLayout) -> None:
        intro = QLabel(
            "Choose one location source. City and Search results are mutually exclusive.",
            tab_widget,
        )
        intro.setWordWrap(True)
        layout.addRow(intro)

        radio_row = QWidget(tab_widget)
        radio_layout = QHBoxLayout(radio_row)
        radio_layout.setContentsMargins(0, 0, 0, 0)
        radio_layout.setSpacing(12)
        self._location_city_radio = QRadioButton("City", tab_widget)
        self._location_place_radio = QRadioButton("Search results", tab_widget)
        self._location_mode_button_group.addButton(self._location_city_radio)
        self._location_mode_button_group.addButton(self._location_place_radio)
        self._location_city_radio.toggled.connect(
            lambda checked: self._set_location_mode("city", checked)
        )
        self._location_place_radio.toggled.connect(
            lambda checked: self._set_location_mode("place", checked)
        )
        radio_layout.addWidget(self._location_city_radio)
        radio_layout.addWidget(self._location_place_radio)
        radio_layout.addStretch(1)
        layout.addRow("Location source", radio_row)

    def _build_time_tab(self, tab_widget: QWidget, layout: QVBoxLayout) -> None:
        source_label = QLabel("Time source", tab_widget)
        source_row = QWidget(tab_widget)
        source_layout = QHBoxLayout(source_row)
        source_layout.setContentsMargins(0, 0, 0, 0)
        source_layout.setSpacing(12)
        self._time_current_radio = QRadioButton("Current time", tab_widget)
        self._time_relative_radio = QRadioButton("Relative time", tab_widget)
        self._time_absolute_radio = QRadioButton("Absolute time", tab_widget)
        self._time_source_button_group.addButton(self._time_current_radio)
        self._time_source_button_group.addButton(self._time_relative_radio)
        self._time_source_button_group.addButton(self._time_absolute_radio)
        self._time_current_radio.toggled.connect(
            lambda checked: self._set_time_source("current", checked)
        )
        self._time_relative_radio.toggled.connect(
            lambda checked: self._set_time_source("relative", checked)
        )
        self._time_absolute_radio.toggled.connect(
            lambda checked: self._set_time_source("absolute", checked)
        )
        source_layout.addWidget(self._time_current_radio)
        source_layout.addWidget(self._time_relative_radio)
        source_layout.addWidget(self._time_absolute_radio)
        source_layout.addStretch(1)
        layout.addWidget(source_label)
        layout.addWidget(source_row)

        shift_container = QWidget(tab_widget)
        shift_layout = QFormLayout(shift_container)
        shift_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
        shift_layout.setFormAlignment(
            Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
        )
        shift_layout.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )

        absolute_container = QWidget(tab_widget)
        absolute_layout = QFormLayout(absolute_container)
        absolute_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
        absolute_layout.setFormAlignment(
            Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
        )
        absolute_layout.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )

        self._time_group_layouts["shift"] = shift_layout
        self._time_group_layouts["absolute"] = absolute_layout
        self._time_group_widgets["shift"] = []
        self._time_group_widgets["absolute"] = []

        layout.addWidget(shift_container)
        layout.addWidget(absolute_container)
        layout.addStretch(1)

    def _location_mode_from_profile(self, profile: dict[str, Any]) -> str:
        city_value = profile.get("city")
        if isinstance(city_value, dict):
            resolver = str(city_value.get("resolver", "")).strip().lower()
            if resolver == "nominatim":
                return "place"
        return "city"

    def _set_location_group_enabled(self, group: str, enabled: bool) -> None:
        for widget in self._location_group_widgets.get(group, []):
            widget.setEnabled(enabled)

    def _apply_location_mode(self, mode: str) -> None:
        if self._location_city_radio is None or self._location_place_radio is None:
            return
        city_checked = mode == "city"
        place_checked = mode == "place"
        city_blocked = self._location_city_radio.blockSignals(True)
        place_blocked = self._location_place_radio.blockSignals(True)
        try:
            self._location_city_radio.setChecked(city_checked)
            self._location_place_radio.setChecked(place_checked)
        finally:
            self._location_city_radio.blockSignals(city_blocked)
            self._location_place_radio.blockSignals(place_blocked)
        self._set_location_group_enabled("city", city_checked)
        self._set_location_group_enabled("place", place_checked)
        if self._city_auto_button is not None:
            self._city_auto_button.setEnabled(city_checked)
        if self._city_search_button is not None:
            self._city_search_button.setEnabled(place_checked)

    def _set_location_mode(self, mode: str, checked: bool) -> None:
        if checked:
            self._apply_location_mode(mode)
            return
        other_mode = "place" if mode == "city" else "city"
        other_radio = (
            self._location_place_radio if other_mode == "place" else self._location_city_radio
        )
        if other_radio is not None and other_radio.isChecked():
            self._apply_location_mode(other_mode)
        else:
            self._apply_location_mode(mode)

    def _city_widget(self) -> QPlainTextEdit | None:
        widget = self._widgets.get("city")
        if isinstance(widget, QPlainTextEdit):
            return widget
        return None

    def _get_city_text(self) -> str:
        widget = self._city_widget()
        if widget is None:
            return ""
        return widget.toPlainText().strip()

    def _set_city_text(self, text: str) -> None:
        widget = self._city_widget()
        if widget is None:
            return
        blocked = widget.blockSignals(True)
        try:
            widget.setPlainText(text)
        finally:
            widget.blockSignals(blocked)

    def _place_payload_from_target(
        self,
        query: str,
        countrycode: str | None,
        language: str,
        target: object,
    ) -> dict[str, Any] | None:
        label = getattr(target, "object_key", "") or getattr(target, "label", "")
        if not isinstance(label, str) or not label.strip():
            return None
        latitude = getattr(target, "latitude_deg", None)
        longitude = getattr(target, "longitude_deg", None)
        if not isinstance(latitude, (int, float)) or not isinstance(longitude, (int, float)):
            return None
        payload: dict[str, Any] = {
            "resolver": "nominatim",
            "query": query,
            "countrycode": countrycode,
            "language": language,
            "display_name": label.strip(),
            "result": {
                "name": label.strip(),
                "lat": float(latitude),
                "lon": float(longitude),
            },
        }
        subtitle = getattr(target, "subtitle", "")
        if isinstance(subtitle, str) and subtitle.strip():
            payload["subtitle"] = subtitle.strip()
        return payload

    def _set_place_selection(self, payload: dict[str, Any] | None) -> None:
        self._place_selected_payload = payload

    def _clear_place_selection(self) -> None:
        self._set_place_selection(None)

    def _build_place_result_payload(
        self,
        query: str,
        countrycode: str | None,
        language: str,
        target: object,
    ) -> dict[str, Any] | None:
        payload = self._place_payload_from_target(query, countrycode, language, target)
        if payload is None:
            return None
        return payload

    def _open_place_search_dialog(self) -> None:
        def _search_callback(search_query: str, countrycode: str | None, language: str):
            candidates = search_place_candidates(
                search_query,
                countrycode=countrycode,
                language=language,
            )
            return build_place_search_jump_targets(candidates)

        dialog = PlaceSearchDialog(
            _search_callback,
            parent=self,
            initial_query="",
            initial_countrycode="",
            initial_language="en",
            show_cli_view_center_checkbox=False,
        )
        dialog.setWindowTitle("Search Place for Startup")
        if dialog.exec() != 1:
            return
        selected = dialog.selected_target()
        if selected is None:
            self._show_validation_error("Select one place candidate before continuing")
            return
        query = dialog.search_query()
        countrycode = dialog.search_countrycode()
        language = dialog.search_language()
        payload = self._build_place_result_payload(query, countrycode, language, selected)
        if payload is None:
            self._show_validation_error("Selected place candidate is invalid")
            return
        self._set_place_selection(payload)
        display_name = str(getattr(selected, "object_key", "") or getattr(selected, "label", "")).strip()
        if display_name:
            self._set_city_text(display_name)

    def _time_source_from_profile(self, profile: dict[str, Any]) -> str:
        hours = float(profile.get("hours", 0.0) or 0.0)
        days = float(profile.get("days", 0.0) or 0.0)
        datetime_text = str(profile.get("datetime", "") or "").strip()
        timezone_text = str(profile.get("timezone", "") or "").strip()
        if datetime_text or timezone_text:
            return "absolute"
        if abs(hours) > 1e-12 or abs(days) > 1e-12:
            return "relative"
        return "current"

    def _set_time_group_enabled(self, group: str, enabled: bool) -> None:
        for widget in self._time_group_widgets.get(group, []):
            widget.setEnabled(enabled)

    def _apply_time_source(self, mode: str) -> None:
        if (
            self._time_current_radio is None
            or self._time_relative_radio is None
            or self._time_absolute_radio is None
        ):
            return
        current_checked = mode == "current"
        relative_checked = mode == "relative"
        absolute_checked = mode == "absolute"
        current_blocked = self._time_current_radio.blockSignals(True)
        relative_blocked = self._time_relative_radio.blockSignals(True)
        absolute_blocked = self._time_absolute_radio.blockSignals(True)
        try:
            self._time_current_radio.setChecked(current_checked)
            self._time_relative_radio.setChecked(relative_checked)
            self._time_absolute_radio.setChecked(absolute_checked)
        finally:
            self._time_current_radio.blockSignals(current_blocked)
            self._time_relative_radio.blockSignals(relative_blocked)
            self._time_absolute_radio.blockSignals(absolute_blocked)
        self._set_time_group_enabled("shift", relative_checked)
        self._set_time_group_enabled("absolute", absolute_checked)

    def _set_time_source(self, mode: str, checked: bool) -> None:
        if checked:
            self._apply_time_source(mode)
            return
        if self._time_current_radio is not None and self._time_current_radio.isChecked():
            self._apply_time_source("current")
            return
        if self._time_relative_radio is not None and self._time_relative_radio.isChecked():
            self._apply_time_source("relative")
            return
        if self._time_absolute_radio is not None and self._time_absolute_radio.isChecked():
            self._apply_time_source("absolute")
        else:
            self._apply_time_source("current")

    def _restore_place_selection(self, profile: dict[str, Any]) -> None:
        city_value = profile.get("city")
        city_widget = self._city_widget()
        if city_widget is None:
            self._set_place_selection(None)
            return
        if not isinstance(city_value, dict):
            self._set_place_selection(None)
            return
        resolver = str(city_value.get("resolver", "")).strip().lower()
        if resolver == "nominatim":
            query = str(city_value.get("query", "")).strip()
            result = city_value.get("result")
            if not isinstance(result, dict):
                self._set_place_selection(None)
                return
            try:
                lat = float(result.get("lat"))
                lon = float(result.get("lon"))
            except (TypeError, ValueError):
                self._set_place_selection(None)
                return
            payload = {
                "resolver": "nominatim",
                "query": query,
                "countrycode": city_value.get("countrycode"),
                "language": str(city_value.get("language", "")).strip() or "en",
                "result": {
                    "name": str(result.get("name", "")).strip(),
                    "lat": lat,
                    "lon": lon,
                },
            }
            if result.get("category") is not None:
                payload["result"]["category"] = result.get("category")
            if result.get("type") is not None:
                payload["result"]["type"] = result.get("type")
            if result.get("importance") is not None:
                payload["result"]["importance"] = result.get("importance")
            subtitle = str(city_value.get("subtitle", "")).strip()
            if subtitle:
                payload["subtitle"] = subtitle
            self._set_place_selection(payload)
            display_name = str(city_value.get("display_name", "")).strip() or str(result.get("name", "")).strip()
            self._set_city_text(display_name)
            return
        if resolver == "auto":
            display_name = str(city_value.get("display_name", "")).strip()
            self._set_city_text(display_name)
            self._set_place_selection(None)
            return
        self._set_place_selection(None)

    def _add_spec(self, spec: _FieldSpec) -> None:
        location_group: str | None = None
        time_group: str | None = None
        if spec.tab == "Time":
            if spec.key in {"hours", "days"}:
                layout = self._time_group_layouts["shift"]
                time_group = "shift"
            elif spec.key in {"datetime", "timezone"}:
                layout = self._time_group_layouts["absolute"]
                time_group = "absolute"
            else:
                raise ValueError(f"Unsupported time field: {spec.key}")
        elif spec.tab == "Overlays":
            section = self._overlay_section_by_key[spec.key]
            layout = self._overlay_layouts[section]
        else:
            layout = self._tab_layouts[spec.tab]
        widget: QWidget
        if spec.kind == "text":
            if spec.key == "city":
                city_edit = QPlainTextEdit(self)
                if spec.placeholder:
                    city_edit.setPlaceholderText(spec.placeholder)
                city_edit.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
                city_edit.setMinimumHeight(40)
                row_widget = QWidget(self)
                row_layout = QVBoxLayout(row_widget)
                row_layout.setContentsMargins(0, 0, 0, 0)
                row_layout.setSpacing(6)
                button_row = QHBoxLayout()
                button_row.setContentsMargins(0, 0, 0, 0)
                button_row.setSpacing(6)
                self._city_auto_button = QPushButton("Auto Search", self)
                self._city_auto_button.clicked.connect(self._start_city_auto)
                button_row.addWidget(self._city_auto_button)
                self._city_search_button = QPushButton("Search ...", self)
                self._city_search_button.clicked.connect(self._open_place_search_dialog)
                button_row.addWidget(self._city_search_button)
                button_row.addStretch(1)
                row_layout.addLayout(button_row)
                row_layout.addWidget(city_edit)
                widget = city_edit
                self._widgets[spec.key] = widget
                self._location_group_widgets["city"].extend([widget, self._city_auto_button])
                self._location_group_widgets["place"].append(self._city_search_button)
                layout.addRow(spec.label, row_widget)
                return
            widget = QLineEdit(self)
            if spec.placeholder:
                widget.setPlaceholderText(spec.placeholder)
        elif spec.kind == "float":
            widget = QDoubleSpinBox(self)
            widget.setDecimals(spec.decimals)
            widget.setRange(spec.minimum, spec.maximum)
            widget.setSingleStep(spec.step)
            widget.setKeyboardTracking(False)
            if spec.key == "view_center_az":
                self._view_center_az_widget = widget
        elif spec.kind == "int":
            widget = QSpinBox(self)
            widget.setRange(int(spec.minimum), int(spec.maximum))
            widget.setSingleStep(int(spec.step))
            widget.setKeyboardTracking(False)
        elif spec.kind == "bool":
            widget = QCheckBox(self)
        elif spec.kind == "choice":
            combo = QComboBox(self)
            combo.addItems(spec.choices)
            widget = combo
        elif spec.kind == "note":
            widget = QLabel(spec.note, self)
            widget.setWordWrap(True)
        else:
            raise ValueError(f"Unsupported field kind: {spec.kind}")
        self._widgets[spec.key] = widget
        if spec.tab == "Location":
            if spec.key == "city":
                location_group = "city"
        elif time_group is not None:
            self._time_group_widgets[time_group].append(widget)
        if location_group is not None:
            self._location_group_widgets[location_group].append(widget)
        if spec.key == "view_center_alt":
            row_widget = QWidget(self)
            row_layout = QVBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(6)
            self._view_center_alt_hint_label = QLabel("Alt value: 0 is horizontal, 90 is zenith.", self)
            self._view_center_alt_hint_label.setWordWrap(True)
            row_layout.addWidget(self._view_center_alt_hint_label)
            row_layout.addWidget(widget)
            layout.addRow(spec.label, row_widget)
            return
        if spec.key == "view_center_az":
            row_widget = QWidget(self)
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(6)
            row_layout.addWidget(widget, 1)
            for label, value in (("N", 0.0), ("E", 90.0), ("S", 180.0), ("W", 270.0)):
                button = QPushButton(label, self)
                button.setFixedWidth(28)
                button.clicked.connect(lambda _checked=False, angle=value: self._set_view_center_az(angle))
                row_layout.addWidget(button)
                self._view_center_az_buttons[label] = button
            layout.addRow(spec.label, row_widget)
            return
        layout.addRow(spec.label, widget)

    def _set_view_center_az(self, value: float) -> None:
        if self._view_center_az_widget is None:
            return
        self._view_center_az_widget.setValue(float(value))

    def _start_city_auto(self) -> None:
        if self._city_auto_button is not None:
            self._city_auto_button.setEnabled(False)
        self._city_auto_request_id += 1
        request_id = self._city_auto_request_id
        submit_gui_work(self._run_city_auto, request_id=request_id)

    def _run_city_auto(self, request_id: int) -> None:
        try:
            resolved = self._auto_location_resolver()
        except Exception as exc:  # pragma: no cover - surfaced through signal delivery
            self.city_auto_finished.emit(request_id, None, f"Auto location failed: {exc}")
            return
        self.city_auto_finished.emit(
            request_id,
            resolved,
            f"Auto location: {resolved.display_name}",
        )

    def _on_city_auto_finished(self, request_id: int, payload: object, status_text: str) -> None:
        if request_id != self._city_auto_request_id:
            return
        if self._city_auto_button is not None:
            self._city_auto_button.setEnabled(True)
        if isinstance(payload, ResolvedLocation):
            self._set_city_text(payload.display_name)
            return
        self._show_error_dialog("Auto location failed", status_text)

    def _restore_from_profile(self, profile: dict[str, Any]) -> None:
        for key, widget in self._widgets.items():
            value = profile.get(key, self._defaults.get(key))
            if key == "city":
                continue
            if isinstance(widget, QLineEdit):
                widget.setText(_as_text(value, key=key))
            elif isinstance(widget, QDoubleSpinBox):
                widget.setValue(
                    _coerce_float_value(key, value, float(widget.value()))
                )
            elif isinstance(widget, QSpinBox):
                widget.setValue(_coerce_int_value(key, value, int(widget.value())))
            elif isinstance(widget, QCheckBox):
                widget.setChecked(bool(value))
            elif isinstance(widget, QComboBox):
                if key in {"show_dso_initial", "show_asterisms_initial", "show_guidelines_initial"}:
                    text = "default" if value is None else ("true" if bool(value) else "false")
                else:
                    text = str(value if value is not None else self._defaults.get(key, ""))
                index = widget.findText(text)
                widget.setCurrentIndex(max(0, index))
        city_value = profile.get("city")
        if isinstance(city_value, dict):
            resolver = str(city_value.get("resolver", "")).strip().lower()
            if resolver == "auto":
                display_name = str(city_value.get("display_name", "")).strip()
                if display_name:
                    self._set_city_text(display_name)
                self._clear_place_selection()
            elif resolver == "nominatim":
                self._restore_place_selection(profile)
            else:
                self._clear_place_selection()
        else:
            self._set_city_text(_as_text(city_value, key="city"))
            self._clear_place_selection()
        self._apply_location_mode(self._location_mode_from_profile(profile))
        self._apply_time_source(self._time_source_from_profile(profile))

    def reset_to_defaults(self) -> None:
        self._base_profile = dict(self._defaults)
        self._restore_from_profile(self._defaults)

    def selected_profile(self) -> dict[str, Any]:
        profile = dict(self._base_profile)
        for key, widget in self._widgets.items():
            if key == "city":
                profile[key] = self._get_city_text()
                continue
            if isinstance(widget, QLineEdit):
                text = widget.text().strip()
                if key == "window_geometry":
                    profile[key] = text or "restore"
                elif key == "cloud_stripe":
                    profile[key] = text
                else:
                    profile[key] = text
            elif isinstance(widget, QDoubleSpinBox):
                profile[key] = float(widget.value())
            elif isinstance(widget, QSpinBox):
                profile[key] = int(widget.value())
            elif isinstance(widget, QCheckBox):
                profile[key] = bool(widget.isChecked())
            elif isinstance(widget, QComboBox):
                text = widget.currentText().strip()
                if key in {"show_dso_initial", "show_asterisms_initial", "show_guidelines_initial"}:
                    profile[key] = _tri_bool_to_value(text)
                else:
                    profile[key] = text
        if self._location_city_radio is not None and self._location_city_radio.isChecked():
            self._clear_place_selection()
        elif self._location_place_radio is not None and self._location_place_radio.isChecked():
            if self._place_selected_payload is not None:
                profile["city"] = self._place_selected_payload
            else:
                profile["city"] = ""
        time_source = "current"
        if self._time_relative_radio is not None and self._time_relative_radio.isChecked():
            time_source = "relative"
        elif self._time_absolute_radio is not None and self._time_absolute_radio.isChecked():
            time_source = "absolute"
        if time_source == "relative":
            profile["datetime"] = None
            profile["timezone"] = None
        elif time_source == "absolute":
            profile["hours"] = 0.0
            profile["days"] = 0.0
        else:
            profile["hours"] = 0.0
            profile["days"] = 0.0
            profile["datetime"] = None
            profile["timezone"] = None
        return profile

    def accept(self) -> None:
        try:
            self._validate_profile(self.selected_profile())
        except ValueError as exc:
            self._show_validation_error(str(exc))
            return
        super().accept()

    def _validate_profile(self, profile: dict[str, Any]) -> None:
        window_geometry = str(profile.get("window_geometry", "")).strip()
        if window_geometry:
            try:
                _parse_window_geometry(window_geometry)
            except Exception as exc:
                raise ValueError("Invalid window geometry") from exc
        cloud_stripe = str(profile.get("cloud_stripe", "")).strip()
        if cloud_stripe:
            try:
                _parse_cloud_stripe(cloud_stripe)
            except Exception as exc:
                raise ValueError("Invalid cloud stripe value") from exc
        if self._location_place_radio is not None and self._location_place_radio.isChecked():
            city_payload = profile.get("city")
            if not isinstance(city_payload, dict) or str(city_payload.get("resolver", "")).strip().lower() != "nominatim":
                raise ValueError("Search and select a place candidate before confirming")
        time_source = "current"
        if self._time_relative_radio is not None and self._time_relative_radio.isChecked():
            time_source = "relative"
        elif self._time_absolute_radio is not None and self._time_absolute_radio.isChecked():
            time_source = "absolute"
        hours = float(profile.get("hours", 0.0) or 0.0)
        days = float(profile.get("days", 0.0) or 0.0)
        datetime_text = str(profile.get("datetime", "") or "").strip()
        has_relative_shift = abs(hours) > 1e-12 or abs(days) > 1e-12
        if time_source == "relative":
            if not has_relative_shift:
                raise ValueError("Relative time requires Hours or Days")
            if datetime_text:
                raise ValueError("Relative shift and absolute time cannot be used together")
        elif time_source == "absolute":
            if not datetime_text:
                raise ValueError("Date/time is required when Absolute time is selected")

    def _show_validation_error(self, message: str) -> None:
        self._show_error_dialog("Invalid input", message)

    def _show_error_dialog(self, title: str, message: str) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        dialog_layout = QVBoxLayout(dialog)
        label = QLabel(message, dialog)
        label.setWordWrap(True)
        dialog_layout.addWidget(label)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok, parent=dialog)
        button_box.accepted.connect(dialog.accept)
        dialog_layout.addWidget(button_box)
        dialog.exec()


def _resolve_auto_location_for_dialog() -> ResolvedLocation:
    return resolve_launch_location(
        "auto",
        load_last_city_func=lambda: None,
        save_last_city_func=lambda _value: None,
    )
