from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QFrame,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ..cli.args import (
    _parse_cloud_stripe,
    _parse_window_geometry,
)
from ..paths import (
    OVERLAY_FONT_SIZE_MAX,
    OVERLAY_FONT_SIZE_MIN,
)
from .launch_profile import default_gui_launch_profile

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
    decimals: int = 2
    minimum: float = 0.0
    maximum: float = 1.0
    step: float = 1.0
    placeholder: str = ""
    note: str = ""


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
    def __init__(
        self,
        profile: dict[str, Any] | None = None,
        *,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("zstarview-gui startup")
        self.setModal(True)
        self.resize(920, 760)

        self._defaults = default_gui_launch_profile()
        self._base_profile = dict(self._defaults)
        if profile:
            self._base_profile.update(profile)

        self._widgets: dict[str, Any] = {}
        self._tab_layouts: dict[str, QFormLayout] = {}

        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(12, 12, 12, 12)

        self._tabs = QTabWidget(self)
        outer_layout.addWidget(self._tabs, 1)

        tab_order = (
            "Location & Time",
            "Stars",
            "Sky",
            "Overlays",
            "General",
            "Search Objects at Startup",
        )
        for tab_name in tab_order:
            tab_widget = QWidget(self)
            tab_layout = QFormLayout(tab_widget)
            tab_layout.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
            tab_layout.setFormAlignment(
                Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
            )
            tab_layout.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
            self._tab_layouts[tab_name] = tab_layout
            scroll_area = QScrollArea(self)
            scroll_area.setWidgetResizable(True)
            scroll_area.setFrameShape(QFrame.Shape.NoFrame)
            scroll_area.setWidget(tab_widget)
            self._tabs.addTab(scroll_area, tab_name)

        self._add_specs()
        self._restore_from_profile(self._base_profile)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self._reset_button = QPushButton("Reset", self)
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
            _FieldSpec("city", "City", "text", "Location & Time"),
            _FieldSpec("place", "Place query", "text", "Location & Time"),
            _FieldSpec(
                "place_countrycode",
                "Place country code",
                "text",
                "Location & Time",
            ),
            _FieldSpec("place_lang", "Place language", "text", "Location & Time"),
            _FieldSpec("hours", "Hours", "float", "Location & Time", minimum=-9999.0, maximum=9999.0, step=0.5),
            _FieldSpec("days", "Days", "float", "Location & Time", minimum=-9999.0, maximum=9999.0, step=0.5),
            _FieldSpec("datetime", "Date/time", "text", "Location & Time"),
            _FieldSpec("timezone", "Timezone", "text", "Location & Time"),
            _FieldSpec("observer_height_m", "Observer height", "float", "Location & Time", minimum=0.0, maximum=10000.0, step=0.1),
            _FieldSpec("use_building_top", "Use building top", "bool", "Location & Time"),
            _FieldSpec("view_center_alt", "View center alt", "float", "Location & Time", minimum=-90.0, maximum=90.0, step=1.0),
            _FieldSpec("view_center_az", "View center az", "float", "Location & Time", minimum=0.0, maximum=360.0, step=1.0),
            _FieldSpec("edge_fov_deg", "Edge FOV", "float", "Location & Time", minimum=0.1, maximum=135.0, step=0.5),
            _FieldSpec("content_fov_deg", "Content FOV", "float", "Location & Time", minimum=90.0, maximum=135.0, step=0.5),
            _FieldSpec("sky_opacity", "Sky opacity", "float", "Sky", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("sky_disc_style", "Sky disc style", "choice", "Sky", choices=("smooth",)),
            _FieldSpec("sky_disc_altaz_rings", "Sky disc rings", "choice", "Sky", choices=("off", "dimalt", "altaz")),
            _FieldSpec("sky_disc_altaz_rings_hover", "Hover rings", "choice", "Sky", choices=("off", "dimalt", "altaz")),
            _FieldSpec("vmag_limit", "Vmag limit", "float", "Stars", minimum=0.0, maximum=20.0, step=0.1),
            _FieldSpec("vmag_brightness_multiplier", "Brightness multiplier", "float", "Stars", minimum=1.0, maximum=3.0, step=0.01),
            _FieldSpec("enlarge_moon", "Enlarge moon", "bool", "Stars"),
            _FieldSpec("bright_bodies", "Bright bodies", "choice", "Stars", choices=("outline", "fill")),
            _FieldSpec("star_base_radius", "Star base radius", "float", "Stars", minimum=0.0, maximum=20.0, step=0.1),
            _FieldSpec("expected_render_width", "Expected render width", "int", "Stars", minimum=1.0, maximum=10000.0, step=1.0),
            _FieldSpec("show_dso_initial", "DSO visibility", "choice", "Stars", choices=("default", "true", "false")),
            _FieldSpec("show_asterisms_initial", "Asterisms visibility", "choice", "Stars", choices=("default", "true", "false")),
            _FieldSpec("show_guidelines_initial", "Guidelines visibility", "choice", "Stars", choices=("default", "true", "false")),
            _FieldSpec("cloud_opacity", "Cloud opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("cloud_stripe", "Cloud stripe", "text", "Overlays"),
            _FieldSpec("cloud_missing_tint_opacity", "Cloud missing tint", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("aircraft_opacity", "Aircraft opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("satellite_opacity", "Satellite opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("terrain_horizon_opacity", "Terrain horizon opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.001),
            _FieldSpec("earth_guide_opacity", "Earth guide opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.001),
            _FieldSpec("night_light_opacity", "Night light opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("ground_tint_opacity", "Ground tint opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("overlay_font_size", "Overlay font size", "float", "Overlays", minimum=float(OVERLAY_FONT_SIZE_MIN), maximum=float(OVERLAY_FONT_SIZE_MAX), step=0.5),
            _FieldSpec("urban_outline_opacity", "Urban outline opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("water_surface_opacity", "Water surface opacity", "float", "Overlays", minimum=0.0, maximum=1.0, step=0.01),
            _FieldSpec("urban_outline_feature_type", "Urban outline mode", "choice", "Overlays", choices=("both", "building")),
            _FieldSpec("urban_outline_radius_km", "Urban radius km", "float", "Overlays", minimum=0.0, maximum=1000.0, step=0.1),
            _FieldSpec("urban_outline_skyscraper_radius_km", "Skyscraper radius km", "float", "Overlays", minimum=0.0, maximum=1000.0, step=0.1),
            _FieldSpec("urban_outline_min_height_m", "Min building height", "float", "Overlays", minimum=0.0, maximum=100000.0, step=0.1),
            _FieldSpec("urban_outline_skyscraper_only", "Skyscraper only", "bool", "Overlays"),
            _FieldSpec("theme", "Theme", "choice", "General", choices=("night", "day", "white", "black", "transparent")),
            _FieldSpec("window_geometry", "Window geometry", "text", "General"),
            _FieldSpec("window_frame", "Window frame", "choice", "General", choices=("frameless", "window")),
            _FieldSpec("observation_info", "Observation info", "choice", "General", choices=("auto", "top", "bottom", "off")),
            _FieldSpec("visibility_boost", "Visibility boost", "float", "General", minimum=1.0, maximum=10.0, step=0.1),
            _FieldSpec("search", "Search query", "text", "Search Objects at Startup"),
            _FieldSpec("search_keep_marker", "Keep marker", "bool", "Search Objects at Startup"),
        )
        for spec in specs:
            self._add_spec(spec)

    def _add_spec(self, spec: _FieldSpec) -> None:
        layout = self._tab_layouts[spec.tab]
        widget: QWidget
        if spec.kind == "text":
            widget = QLineEdit(self)
            if spec.placeholder:
                widget.setPlaceholderText(spec.placeholder)
        elif spec.kind == "float":
            widget = QDoubleSpinBox(self)
            widget.setDecimals(spec.decimals)
            widget.setRange(spec.minimum, spec.maximum)
            widget.setSingleStep(spec.step)
            widget.setKeyboardTracking(False)
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
        layout.addRow(spec.label, widget)

    def _restore_from_profile(self, profile: dict[str, Any]) -> None:
        for key, widget in self._widgets.items():
            value = profile.get(key, self._defaults.get(key))
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

    def reset_to_defaults(self) -> None:
        self._base_profile = dict(self._defaults)
        self._restore_from_profile(self._defaults)

    def selected_profile(self) -> dict[str, Any]:
        profile = dict(self._base_profile)
        for key, widget in self._widgets.items():
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

    def _show_validation_error(self, message: str) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("Invalid input")
        dialog_layout = QVBoxLayout(dialog)
        label = QLabel(message, dialog)
        label.setWordWrap(True)
        dialog_layout.addWidget(label)
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok, parent=dialog)
        button_box.accepted.connect(dialog.accept)
        dialog_layout.addWidget(button_box)
        dialog.exec()
