# レンダリング・パイプライン

この文書は、`zstarview`、`zstarview-gui`、`zstarview-export-image` で共通する描画パイプラインをまとめる。
主要なデータモデルは [data-model.md](data-model.md) を参照する。

## 1. 役割

- GUI と export-image が同じ入力から同じ scene を作れるよう、共通パイプラインを使う。
- ベースシーン描画と HUD / 一時オーバーレイ描画を分ける。
- ラベル、ガイド、オーバーレイは、ウィジェットへ直結させず、合成可能な部品として扱う。
- renderer 固有の ad hoc な判断より、前処理済みの scene data と state を優先する。

## 2. 合成の順序

描画は、概ね次の層を順に重ねる。

1. 背景
2. sky disc
3. 雲と欠損ティント
4. sky guides
5. 地形地平線と地面ティント
6. Earth guide
7. 星と天体
8. 都市アウトライン
9. 水面、航空機、人工衛星、夜間光、台風・サイクロン
10. ラベルと HUD

この順序は厳密な固定ではなく、レイヤーの種類ごとの責務を表す。
たとえば、雲は sky disc の上に乗るが、地平線下のティントや HUD はさらに上に置く。

## 3. ベース描画と HUD

- ベース描画は、scene の視覚的な土台を作る。
- HUD は、hover、jump highlight、static observation overlay、status line のような一時的または情報密度の高い表示を担当する。
- ベース描画と HUD を分けることで、視点変更時に base frame を再利用しやすくなる。
- `zstarview-export-image` は、既定では HUD を入れずにベース描画中心で出力してよい。

### 3.1 star-focus mode

- `star-focus` mode は、通常の HUD を抑えたまま、星と主要なガイドだけを見せるための GUI 専用の表示モードとする。
- このモードは `viewport_interaction_mode` と独立させる。
- `viewport_interaction_mode` は矢印キー操作や resize 時の一時軽量化を担い、hover や HUD を含めた通常の一時制御とは別責務にする。
- `star-focus` mode は hover 判定を維持し、hover ラベルは残してよい。
- `star-focus` mode では、status line と observation HUD を非表示にしてよい。
- クリックで入る起動条件、再クリックによる解除、ドラッグ開始時の解除は GUI 層で状態を更新し、renderer には最終的な active 状態だけ渡してよい。

## 4. オーバーレイの責務

### 4.1 天体とガイド

- 星、太陽、月、惑星は主要天体として描画する。
- Sky guides は、地平線、天の赤道、黄道、方位ラベル、天頂マーカー、天の極マーカーを含む。
- never-rises 領域は、Sky guides の一部として扱ってよい。

### 4.2 雲と地平線

- 雲は sky disc の上に合成し、欠損領域は別ティントで表現してよい。
- 地形地平線は空と地面の境界を決める。
- Earth guide は地平線下の方向案内として、地形地平線とは独立した補助レイヤーとする。

### 4.3 地表系レイヤー

- 都市アウトラインは、観測地点周辺の建物輪郭を示す。
- 水面レイヤーは、地形地平線や地面ティントと合わせて地表の読み取りを助ける。

### 4.4 移動体と補助情報

- 航空機は予測軌跡として描画する。
- 人工衛星は小さなマーカーとして描画する。
- 夜間光と台風・サイクロンは、補助的な地理・気象情報として重ねる。
- 夜間光のうち、主稜線に沿う sky glow は `src/zstarview/render/ridge_glow.py` の独立した補助レイヤーとして扱う。街灯系の glow と色推定は分けてよい。

## 5. ラベルの扱い

- ラベルは、星、惑星、DSO、アステリズム、衛星、航空機、検索結果など複数の候補源から集約する。
- まず候補を正規化し、その後に配置順を決める。
- 重なりやすいラベルは、局所グループを作ってから配置してよい。
- グループ内では、対象に近い候補、または上側の候補から先に配置してよい。
- 既に置かれたラベルの予約矩形を共有し、後続のラベルほど回避するようにしてよい。
- hover ラベルと persistent marker ラベルは、短い jump highlight とは別に扱ってよい。
- export-image では、必要なら marker/label を 1 件だけ持ち込むが、永続状態は残さなくてよい。

## 6. オーバーレイと HUD の境界

- static observation overlay は HUD 側で描画する。
- status line も HUD 側で描画する。
- hover の一時表示は、ベース描画ではなく HUD 側で足す。
- ベースフレームのキャッシュから mouse position、hover 対象、jump highlight、status 文言は除外してよい。
- `star-focus` mode のときは、HUD 側のうち location / time / status / observation block を抑止しつつ、hover に必要な最小限の描画とラベル候補の解決は残してよい。

## 7. 外部依存

- GUI toolkit と描画プリミティブ
- 天文・暦関連ライブラリ
- 数値処理
- ラスタ・地理空間処理
- 衛星、地形、建物、航空機、夜間光の取得ヘルパ

## 8. 設計意図

renderer は、正規化済みの scene data を受け取り、入口ごとの分岐を最小限にして最終画像を生成する。
それにより、3 つのアプリの見え方を揃えつつ、起動方法や永続化だけを入口ごとに変えられる。
