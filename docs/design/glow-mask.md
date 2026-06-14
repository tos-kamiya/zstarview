# GlowMask 設計

## 1. 目的

GlowMask は、夜間光と ridge glow の「にじみ」だけを別経路で保持する中間表現である。
狙いは次のとおり。

- 帯ポリゴンの積み上げで発生する banding を避ける。
- 光の強さを float のまま最後まで保持する。
- 色は最後に一度だけ tint し、alpha だけで減衰を表す。
- GUI と export-image で同じ合成結果を得る。

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
3. 低解像度キャンバスに、night-light glow と ridge glow を同じ面へ描く。
   - 帯の内部は単色ではなく、下端から上端へ alpha が連続的に変わるグラデーションで埋める。
4. 低解像度 RGBA から alpha 面だけを抜き出す。
5. 小さな平滑化を数パスかける。
6. `GlowMask` として返す。

この経路では、街灯系 glow と ridge glow は「同じ発光場の別寄与」として扱う。
個別の帯や層を保持するのではなく、最終的な alpha 勾配へ畳み込む。

## 4. ラスタ化の原則

- 低解像度化は `scale` で一律に行う。
- クリップは地平線外へはみ出さないよう、`content_fov_deg / edge_fov_deg` に応じた円形領域で行う。
- 既存の night-light / ridge glow の描画関数を再利用し、GlowMask 用の特別な幾何計算を増やしすぎない。
- `fast_mode` では GlowMask を完全にスキップする。

GlowMask が受け取る描画は、見た目の完成図ではなく「発光の密度場」である。
そのため、細かい帯の境界よりも、面としての連続性を優先する。

## 5. ぼかし

GlowMask の alpha は、低解像度のまま軽くぼかす。

- ぼかしは小さなカーネルを複数回適用する。
- 目的は輪郭の硬さを取ることであり、大きな拡散ではない。
- 強すぎるぼかしは暗く見えるので、パス数と mask scale の両方で調整する。

現在の既定値は次のとおり。

- `GLOW_MASK_SCALE = 0.25`
- `GLOW_MASK_BLUR_PASSES = 2`

## 6. tint と合成

GlowMask は、復元時に固定色へ tint してから RGBA に戻す。

- tint 色は `GLOW_MASK_TINT_RGB` を使う。
- tint の前に HSV の value を最大化し、alpha だけを fade の主体にする。
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
- `ridge_glow_opacity`
- `night_light_sun_alt_deg`
- `fast_mode`
- GlowMask の `scale`
- blur のパラメータ
- tint 色

これにより、見た目に影響する入力が変わったときだけ再計算できる。

## 8. 失敗時の扱い

- 入力サイズが不正なら `None` を返す。
- night-light profile が空、または opacity が両方 0 なら `None` を返す。
- 低解像度マスクに有効な alpha が残らなければ `None` を返す。

つまり GlowMask は「常に何かを描く」前提ではなく、寄与が無いときは無視してよい。

## 9. テスト方針

最低限、次を検証する。

- alpha のぼかしでピークが平坦化すること
- tint が brightness-maximized な base color を使うこと
- 低解像度の描画が geometry に追従すること
- `fast_mode` では GlowMask を経由しないこと

## 10. 今後の拡張

GlowMask を中心に寄せると、今後は次の拡張がしやすい。

- 夜間光と ridge glow の寄与比を分離して調整する
- blur の強さを視野や解像度に応じて可変にする
- 将来的に複数の glow 種別を同じ alpha field に追加する
- 低解像度のまま、より連続的な発光表現へ移行する
