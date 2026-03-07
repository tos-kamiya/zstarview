from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional


@dataclass(frozen=True)
class Asterism:
    key: str
    name: str
    season: str
    edges: tuple[tuple[str, str], ...]

    def segments(self) -> Iterable[tuple[str, str]]:
        return self.edges

    @property
    def source_ids(self) -> tuple[str, ...]:
        ordered: list[str] = []
        seen: set[str] = set()
        for source_a, source_b in self.edges:
            if source_a not in seen:
                seen.add(source_a)
                ordered.append(source_a)
            if source_b not in seen:
                seen.add(source_b)
                ordered.append(source_b)
        return tuple(ordered)


def _hip(hip_id: int) -> str:
    return f"HIP{int(hip_id)}"


def _path_edges(*source_ids: str) -> tuple[tuple[str, str], ...]:
    if len(source_ids) < 2:
        return ()
    return tuple((source_ids[i], source_ids[i + 1]) for i in range(len(source_ids) - 1))


ASTERISMS: tuple[Asterism, ...] = (
    # Winter
    Asterism("winter_triangle", "Winter Triangle", "winter", _path_edges(_hip(32349), _hip(37279), _hip(27989), _hip(32349))),
    Asterism("orions_belt", "Orion's Belt", "winter", _path_edges(_hip(26727), _hip(26311), _hip(25930))),
    Asterism(
        "winter_hexagon",
        "Winter Hexagon",
        "winter",
        _path_edges(_hip(24608), _hip(21421), _hip(24436), _hip(32349), _hip(37279), _hip(37826), _hip(24608)),
    ),
    Asterism(
        "southern_cross",
        "Southern Cross",
        "winter",
        ((_hip(61084), _hip(60718)), (_hip(59747), _hip(62434))),
    ),
    Asterism("southern_pointers", "Southern Pointers", "winter", ((_hip(71683), _hip(68702)),)),
    Asterism("diamond_cross", "Diamond Cross", "winter", ((_hip(45238), _hip(52419)), (_hip(48002), _hip(50099)))),
    Asterism("false_cross", "False Cross", "winter", ((_hip(42913), _hip(45556)), (_hip(45941), _hip(41037)))),
    # Spring
    Asterism(
        "big_dipper",
        "Big Dipper",
        "spring",
        _path_edges(_hip(54061), _hip(53910), _hip(58001), _hip(59774), _hip(62956), _hip(65378), _hip(67301)),
    ),
    Asterism(
        "little_dipper",
        "Little Dipper",
        "spring",
        _path_edges(_hip(11767), _hip(85822), _hip(82080), _hip(77055), _hip(72607), _hip(75097), _hip(79822), _hip(77055)),
    ),
    Asterism("spring_triangle", "Spring Triangle", "spring", _path_edges(_hip(69673), _hip(65474), _hip(57632), _hip(69673))),
    Asterism("arc_to_arcturus", "Arc to Arcturus", "spring", _path_edges(_hip(67301), _hip(69673), _hip(65474))),
    Asterism("leo_sickle", "Leo Sickle", "spring", _path_edges(_hip(49669), _hip(49583), _hip(50583), _hip(50335), _hip(48455), _hip(47908))),
    Asterism("southern_triangle", "Southern Triangle", "spring", _path_edges(_hip(74946), _hip(82273), _hip(77952), _hip(74946))),
    # Summer
    Asterism("summer_triangle", "Summer Triangle", "summer", _path_edges(_hip(91262), _hip(102098), _hip(97649), _hip(91262))),
    Asterism(
        "northern_cross",
        "Northern Cross",
        "summer",
        _path_edges(_hip(102098), _hip(100453), _hip(95947), _hip(100453), _hip(97165), _hip(100453), _hip(102488)),
    ),
    Asterism(
        "teapot",
        "Teapot",
        "summer",
        (
            (_hip(93864), _hip(92855)),
            (_hip(93864), _hip(93506)),
            (_hip(93506), _hip(92104)),
            (_hip(92104), _hip(92855)),
            (_hip(92104), _hip(90496)),
            (_hip(92855), _hip(89931)),
            (_hip(90496), _hip(89931)),
            (_hip(89931), _hip(88635)),
            (_hip(88635), _hip(90185)),
            (_hip(90185), _hip(89931)),
            (_hip(93506), _hip(90185)),
        ),
    ),
    Asterism("keystone", "Keystone", "summer", _path_edges(_hip(84380), _hip(81833), _hip(81693), _hip(83207), _hip(84380))),
    Asterism(
        "coathanger",
        "Coathanger",
        "summer",
        _path_edges(_hip(94703), _hip(95498), _hip(96275), _hip(96757), _hip(97365), _hip(96275), _hip(96837), _hip(96516)),
    ),
    Asterism(
        "scorpion_hook",
        "Scorpion Hook",
        "summer",
        _path_edges(_hip(78820), _hip(78401), _hip(80763), _hip(80112), _hip(82396), _hip(86228), _hip(85927)),
    ),
    # Autumn
    Asterism(
        "great_square_of_pegasus",
        "Great Square of Pegasus",
        "autumn",
        _path_edges(_hip(113963), _hip(113881), _hip(677), _hip(1067), _hip(113963)),
    ),
    Asterism(
        "circlet_of_pisces",
        "Circlet of Pisces",
        "autumn",
        _path_edges(_hip(115830), _hip(117375), _hip(114971), _hip(115738), _hip(116928), _hip(117245), _hip(116771), _hip(115830)),
    ),
    Asterism("water_jar_of_aquarius", "Water Jar of Aquarius", "autumn", _path_edges(_hip(106278), _hip(109074), _hip(110395), _hip(113136), _hip(110003), _hip(109074))),
    Asterism("andromeda_chain", "Andromeda Chain", "autumn", _path_edges(_hip(677), _hip(5447), _hip(9640))),
    Asterism("autumn_triangle", "Autumn Triangle", "autumn", _path_edges(_hip(113963), _hip(9884), _hip(113368), _hip(113963))),
    Asterism("cassiopeia_w", "Cassiopeia W", "autumn", _path_edges(_hip(746), _hip(3179), _hip(4427), _hip(6686), _hip(8886))),
    Asterism(
        "house_of_cepheus",
        "House of Cepheus",
        "autumn",
        _path_edges(_hip(116727), _hip(112724), _hip(106032), _hip(102422), _hip(105199), _hip(109492), _hip(116727)),
    ),
    Asterism("jobs_coffin", "Job's Coffin", "autumn", _path_edges(_hip(101769), _hip(102281), _hip(102532), _hip(101958), _hip(101769))),
)

ASTERISM_KEYS_BY_SOURCE_ID: dict[str, tuple[str, ...]] = {}
for _asterism in ASTERISMS:
    for _source_id in _asterism.source_ids:
        if _source_id not in ASTERISM_KEYS_BY_SOURCE_ID:
            ASTERISM_KEYS_BY_SOURCE_ID[_source_id] = (_asterism.key,)
        else:
            if _asterism.key not in ASTERISM_KEYS_BY_SOURCE_ID[_source_id]:
                ASTERISM_KEYS_BY_SOURCE_ID[_source_id] = ASTERISM_KEYS_BY_SOURCE_ID[_source_id] + (_asterism.key,)

ASTERISM_REQUIRED_SOURCE_IDS: frozenset[str] = frozenset(ASTERISM_KEYS_BY_SOURCE_ID)

ASTERISM_BY_KEY: dict[str, Asterism] = {asterism.key: asterism for asterism in ASTERISMS}


def pick_rotating_asterism(source_id: str, second_slot: int) -> Optional[Asterism]:
    keys = ASTERISM_KEYS_BY_SOURCE_ID.get(str(source_id).strip(), ())
    if not keys:
        return None
    idx = int(second_slot) % len(keys)
    return ASTERISM_BY_KEY.get(keys[idx])
