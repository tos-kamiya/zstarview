• いいえ、現状の ridge glow は night light と同様の alt/az 2D 積算にはなっていませ
  ん。

  主な差分:

  - night light は alpha_grid として alt=-90..90 の 2D field を作っています。
    src/zstarview/night_lights.py:994

  - ridge glow は band_profiles の az-only 強度として作られています。
    src/zstarview/night_lights.py:1025

  - ridge glow は _accumulate_local_glow_strengths() で
    target_altitudes_deg=horizon_alt_values にだけ積算しており、alt 全範囲の field
    には積算していません。
    src/zstarview/night_lights.py:1030

  - ridge glow の source 強度は「一定値のみ」ではなく、sample_matrix +
    NIGHT_LIGHTS_RIDGE_GLOW_SAMPLE_FLOOR です。つまり夜光タイル値も混ざっています。
    src/zstarview/night_lights.py:949

  - ridge glow のマスクは主稜線のみです。ここは意図どおりです。
    src/zstarview/night_lights.py:1038

  - 描画側では ridge glow は band_profiles を平均した az-only brightness に、主稜線
    からの縦方向 falloff をかけています。
    src/zstarview/gui/composite.py:380
    src/zstarview/gui/composite.py:408

  - さらに ridge glow は最終式で brightness、つまり night light 側 brightness にも
    掛けられています。これにより ridge glow は night light field に依存しています。
    src/zstarview/gui/composite.py:419

  結論: ridge glow は「DEM altitude と距離 band を使う」部分は近いですが、night
  light と同じ 2D alt/az field 積算ではありません。
  また、「夜光タイルを使わず一定値」でもなく、現在は夜光タイル値に floor を足してい
  ます。
