# zstarview 設計書 - legacy archive

最終更新: 2026-06-02

この文書は、`docs/design/` に分割する前の単一ファイル版の記録である。
現在の入口は `docs/design.md`、詳細は `docs/design/*.md` を参照する。
これはアーカイブであり、現行設計の正本としては読まなくてよい。

## 1. この文書の位置づけ

この文書は、`zstarview` の内部設計をまとめたものである。  
アーキテクチャ、モジュール責務、主要データ構造、処理フロー、スレッド分離、外部依存の境界を扱う。  
利用者向けの機能説明は `docs/specification.md`、内部判断の履歴は `dev-notes/session-YYYY-MM-DD.md`、過去の実装履歴アーカイブは `docs/implementation-archive.md` を参照する。

## 2. 設計方針

- UI と計算処理を分離する。
- 長時間処理は UI スレッドから切り離す。
- 外部 I/O は局所化し、失敗を UI 側で扱いやすい状態へ正規化する。
- 描画入力はできるだけ前処理済みデータとして渡し、描画側を単純化する。
- CLI 起動時の選択肢と GUI 実行中の状態を同じモデルで扱う。
- 雲、地形地平線、星図更新は相互に独立した補助パイプラインとして扱う。
- 更新頻度の低い補助レイヤーは、連続アニメーションよりも責務分離と UI 安全性を優先する。
- 視点変更の描画は fast-mode と normal-mode を分け、fast-mode は最新 1 件だけ保持する軽い再投影として扱う。
- 定期更新は「前回のダウンロードや計算が終わった時刻」基準で次回期限を決め、worker busy 中は tick を溜めずに次回の idle で 1 件だけ進める。
- ターミナル、console、log、CLI help、例外文言、subprocess の stdout/stderr に出る可能性がある文字列は、Windows 文字コード事故を避けるため ASCII-only を原則とする。
- 非 ASCII 文字を判定ロジックで扱う必要がある場合は、ソースコード中に文字を直書きせず `"\u2019"` のような Unicode escape を優先する。
- UI 専用文字列は、サポート環境で表示確認済みなら非 ASCII を許容する。

## 3. 全体アーキテクチャ

システムは大きく次の層に分かれる。

- 起動・設定層
  - CLI 解析
  - 設定読込
  - 地点解決
  - 初期データロード
- ドメイン計算層
  - 天体位置計算
  - 星カタログ前処理
  - アステリズム、DSO、地形地平線の補助計算
- 描画層
  - 星空ディスク
  - 恒星、惑星、ラベル、補助線描画
  - 雲と地面ティントの合成
  - 地平線下の地球裏面ガイド合成
- UI 層
  - ウィンドウ管理
  - 入力イベント
  - バックグラウンド更新の起動と反映
- 外部データ連携層
  - 衛星クラウドデータ取得
  - DEM 取得
  - 水面レイヤー取得
  - Overture 建物データ取得
  - 夜間光 GeoTIFF 取得
  - 台風・サイクロンの公開 ArcGIS FeatureServer 取得
  - キャッシュ管理

## 4. モジュール構成

### 4.1 起動・設定

- `src/zstarview/gui/viewer.py`
  - GUI アプリケーションの主エントリポイント
  - 起動シーケンスの組み立て
  - `zstarview-gui` の起動前設定ダイアログを起動し、確定値を GUI 初期状態へ反映する
  - startup overlay は初期 sky データと terrain DEM の両方が揃うまで維持し、ready/fail の両経路で初回表示の切り替えを行う
- `src/zstarview/gui/startup_dialog.py`
  - GUI 起動前設定ダイアログを定義する
  - 前回起動値の編集、Reset 操作、確定/取消の結果を返す
  - README の GUI 対応 CLI グループ見出しに対応するタブを持ち、`Location`、`View`、`Time`、`Stars`、`Overlays`、`General` の順に分けて起動前の入力を見通しよく整理する
  - `Overlays` は `Sky`、`Clouds`、`Tropical Cyclone`、`Aircraft and Satellites`、`Ground and Guides`、`Urban Outline` の折りたたみ可能なサブグループへ分けてよい
  - `overlay_font_size` は `General` 側で扱い、全体の表示密度を調整する設定として保持してよい
  - `City` 行には複数行表示の入力欄を置き、その上に `Auto Search` と `Search ...` のボタン行を配置してよい
  - `Auto Search` ボタンは現在地自動取得の結果を `City` 欄へ反映してよい
  - `Location` タブでは `City` と `Search results` を排他的なモードとして扱い、`Search ...` ボタンから専用の place search dialog を開いてよい
  - place search dialog では、自由入力をそのまま Nominatim に投げるのではなく、明示的な検索実行の後に候補一覧から 1 件を選ぶフローとして扱ってよい
  - `Place country code` は空欄なら未指定、`Place language` は空欄なら `en` として正規化し、`None` を Nominatim 呼び出しへ流さない
  - 起動前ダイアログの place 検索は、既存の `place_search_dialog.py` の候補選択フローを再利用してよく、トップ候補の自動採用は避けてよい
