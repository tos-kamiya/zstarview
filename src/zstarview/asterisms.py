from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional


@dataclass(frozen=True)
class Asterism:
    key: str
    name: str
    season: str
    path: tuple[str, ...]

    def segments(self) -> Iterable[tuple[str, str]]:
        if len(self.path) < 2:
            return ()
        return ((self.path[i], self.path[i + 1]) for i in range(len(self.path) - 1))


def _hip(hip_id: int) -> str:
    return f"HIP{int(hip_id)}"


ASTERISMS: tuple[Asterism, ...] = (
    # Winter
    Asterism("winter_triangle", "Winter Triangle", "winter", (_hip(32349), _hip(37279), _hip(27989), _hip(32349))),
    Asterism("orions_belt", "Orion's Belt", "winter", (_hip(26727), _hip(26311), _hip(25930))),
    Asterism("orions_sword", "Orion's Sword", "winter", (_hip(26311), _hip(26241), _hip(27366))),
    Asterism(
        "winter_hexagon",
        "Winter Hexagon",
        "winter",
        (_hip(24608), _hip(21421), _hip(24436), _hip(32349), _hip(37279), _hip(37826), _hip(24608)),
    ),
    Asterism("hyades_v", "Hyades V", "winter", (_hip(20205), _hip(21421), _hip(20889), _hip(20455), _hip(20205))),
    # Spring
    Asterism(
        "big_dipper",
        "Big Dipper",
        "spring",
        (_hip(54061), _hip(53910), _hip(58001), _hip(59774), _hip(62956), _hip(65378), _hip(67301)),
    ),
    Asterism(
        "little_dipper",
        "Little Dipper",
        "spring",
        (_hip(11767), _hip(85822), _hip(82080), _hip(77055), _hip(74793), _hip(75097), _hip(72607), _hip(77055)),
    ),
    Asterism("spring_triangle", "Spring Triangle", "spring", (_hip(69673), _hip(65474), _hip(57632), _hip(69673))),
    Asterism("arc_to_arcturus", "Arc to Arcturus", "spring", (_hip(67301), _hip(65378), _hip(62956), _hip(69673))),
    Asterism("leo_sickle", "Leo Sickle", "spring", (_hip(49669), _hip(50583), _hip(50335), _hip(48455), _hip(47908), _hip(49669))),
    # Summer
    Asterism("summer_triangle", "Summer Triangle", "summer", (_hip(91262), _hip(102098), _hip(97649), _hip(91262))),
    Asterism(
        "northern_cross",
        "Northern Cross",
        "summer",
        (_hip(102098), _hip(100453), _hip(95947), _hip(100453), _hip(97165), _hip(100453), _hip(102488)),
    ),
    Asterism(
        "teapot",
        "Teapot",
        "summer",
        (_hip(88635), _hip(90185), _hip(93506), _hip(92855), _hip(92041), _hip(90496), _hip(89931), _hip(88635)),
    ),
    Asterism("keystone", "Keystone", "summer", (_hip(83207), _hip(81693), _hip(84379), _hip(86974), _hip(83207))),
    Asterism(
        "coathanger",
        "Coathanger",
        "summer",
        (_hip(94703), _hip(95498), _hip(96275), _hip(96757), _hip(97365), _hip(96275), _hip(96837), _hip(96516)),
    ),
    # Autumn
    Asterism(
        "great_square_of_pegasus",
        "Great Square of Pegasus",
        "autumn",
        (_hip(113963), _hip(113881), _hip(677), _hip(1067), _hip(113963)),
    ),
    Asterism("circlet_of_pisces", "Circlet of Pisces", "autumn", (_hip(7097), _hip(8198), _hip(9487), _hip(113889), _hip(109427), _hip(7097))),
    Asterism("water_jar_of_aquarius", "Water Jar of Aquarius", "autumn", (_hip(106278), _hip(109074), _hip(110395), _hip(113136), _hip(110003), _hip(109074))),
    Asterism("andromeda_chain", "Andromeda Chain", "autumn", (_hip(677), _hip(5447), _hip(9640))),
    Asterism("autumn_triangle", "Autumn Triangle", "autumn", (_hip(113963), _hip(9884), _hip(113368), _hip(113963))),
)

ASTERISM_KEYS_BY_SOURCE_ID: dict[str, tuple[str, ...]] = {}
for _asterism in ASTERISMS:
    for _source_id in _asterism.path:
        if _source_id not in ASTERISM_KEYS_BY_SOURCE_ID:
            ASTERISM_KEYS_BY_SOURCE_ID[_source_id] = (_asterism.key,)
        else:
            if _asterism.key not in ASTERISM_KEYS_BY_SOURCE_ID[_source_id]:
                ASTERISM_KEYS_BY_SOURCE_ID[_source_id] = ASTERISM_KEYS_BY_SOURCE_ID[_source_id] + (_asterism.key,)

ASTERISM_BY_KEY: dict[str, Asterism] = {asterism.key: asterism for asterism in ASTERISMS}


def pick_rotating_asterism(source_id: str, second_slot: int) -> Optional[Asterism]:
    keys = ASTERISM_KEYS_BY_SOURCE_ID.get(str(source_id).strip(), ())
    if not keys:
        return None
    idx = int(second_slot) % len(keys)
    return ASTERISM_BY_KEY.get(keys[idx])
