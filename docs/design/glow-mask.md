# GlowMask 設計

## 1. 目的

GlowMask は、夜間光の「にじみ」を別経路で保持する中間表現である。
狙いは次のとおり。

- 帯ポリゴンの積み上げで発生する banding を避ける。
- 光の強さを float のまま最後まで保持する。
- 色は最後に一度だけ tint し、alpha だけで減衰を表す。
- GUI と export-image で同じ合成結果を得る。

街灯系 glow は、低解像度の画面グリッドを逆投影して `alt` / `az` を求め、夜間光プロフィールの azimuth サンプルを補間して alpha 面へ落とす。レイの高さは固定で、地平線から上だけを連続的に減衰させる。ridge glow の簡易シミュレーションでは、profile の最終値を後から持ち上げるのではなく、可視サンプルごとに小さな floor を足してから積算する。

## 2. データモデル

GlowMask は alpha-only の低解像度マスクとして扱う。

- `alpha: np.ndarray`
  - `float32`
  - 値域は `[0.0, 1.0]`
  - 2 次元配列で、mask そのものは色を持たない
- `scale: float`
  - 元の描画解像度に対する縮小率
  - 現行の既定値は `0.25`

GlowMask は scene/state の永続データには入れず、`SkyCompositorCache` の再構築用中間画像として扱う。

## 3. 生成フロー

GlowMask は、最終画像の幅・高さと `ScreenGeometry` から作る。

1. 元画像の `scale` 倍で低解像度の `QImage` を用意する。
2. 低解像度の geometry を作る。
3. 低解像度キャンバスに、night-light glow を描く。
   - 街灯系 glow は、逆投影した `alt` / `az` から night-light profile を補間し、固定高さの連続 alpha field として積み上げる。
   - ridge glow は、街灯よりかなり低い高さの柱として扱い、可視サンプルごとの floor と小さな色混合で horizon 近傍に寄せる。
   - 街灯系 glow のマスクは副稜線レイヤーの「その azimuth までの累積最大 alt」を使う。
   - `terrain_profile_altaz` がない場合は profile の horizon サンプルを使う。
   - 地平線は hard mask ではなく、近傍で滑らかに落ちる soft falloff として扱う。
4. 低解像度 RGBA から alpha 面だけを抜き出す。
5. `GlowMask` として返す。

この経路では、街灯系 glow を単一の発光場として扱う。
個別の帯や層を保持するのではなく、最終的な alpha 勾配へ畳み込む。
さらに、RGBA 化の直前に screen-fixed の deterministic noise を alpha へ掛け、完全に滑らかな見え方を少し崩す。
このノイズはフレームごとに変わらないため、画面に対して安定したザラつきとして見える。

## 4. ラスタ化の原則

- 低解像度化は `scale` で一律に行う。
- クリップは地平線外へはみ出さないよう、`content_fov_deg / edge_fov_deg` に応じた円形領域で行う。
- night-light は逆投影した画面座標上で alpha field を直接作る。
- `fast_mode` では GlowMask を完全にスキップする。

GlowMask が受け取る描画は、見た目の完成図ではなく「発光の密度場」である。
そのため、細かい帯の境界よりも、面としての連続性を優先する。
ただし最終的な見え方は完全な均一面にせず、低解像度の alpha に軽い空間ノイズを乗せて粒立ちを残す。

街灯系 glow の試作では、密度場は「逆投影した画面上の各点が持つ ray-sampled 強度」として扱ってよい。必要な連続性は、azimuth 補間と高さ方向の減衰で確保してよい。ridge glow の擬似分は、可視サンプルごとの小さな floor と、より低い高さの減衰で足し込んでよい。

## 5. tint と合成

GlowMask は、復元時に固定色へ tint してから RGBA に戻す。

- tint 色は `GLOW_MASK_TINT_RGB` を使う。
- tint の前に HSV の value を最大化し、alpha だけを fade の主体にする。
- RGBA 化の直前に、screen-fixed の noise field を alpha へ乗算して軽いザラつきを加える。
- ノイズは低解像度マスク座標から決定論的に生成し、フレーム間で変化しないようにする。
- 出力は premultiplied RGBA の `QImage` にする。
- 最終合成は `CompositionMode_Plus` で行い、他のベース描画を壊さずに重ねる。

この設計では、色味は固定しつつ、発光強度だけを alpha で調整できる。

## 7. キャッシュ

GlowMask は再生成コストを抑えるため、composite キャッシュの一部として扱う。

キャッシュキーには少なくとも次を含める。

- 画面サイズ
- geometry
- view center
- night-light の profile 内容
- `night_light_opacity`
- `night_light_sun_alt_deg`
- `fast_mode`
- GlowMask の `scale`
- tint 色

これにより、見た目に影響する入力が変わったときだけ再計算できる。

## 7. 失敗時の扱い

- 入力サイズが不正なら `None` を返す。
- night-light profile が空、または opacity が両方 0 なら `None` を返す。
- 低解像度マスクに有効な alpha が残らなければ `None` を返す。

つまり GlowMask は「常に何かを描く」前提ではなく、寄与が無いときは無視してよい。

## 8. テスト方針

最低限、次を検証する。

- tint が brightness-maximized な base color を使うこと
- ノイズ付きの alpha 変調がフレーム間で安定していること
- 低解像度の描画が geometry に追従すること
- `fast_mode` では GlowMask を経由しないこと
- night-light の alpha が azimuth 補間と高さ減衰で変化すること

## 9. 今後の拡張

GlowMask を中心に寄せると、今後は次の拡張がしやすい。

- 将来的に複数の glow 種別を同じ alpha field に追加する
- 低解像度のまま、より連続的な発光表現へ移行する
