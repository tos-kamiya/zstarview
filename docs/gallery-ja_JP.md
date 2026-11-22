# zstarview スクリーンショットギャラリー

これらのスクリーンショットは以前のバージョンの zstarview で撮影したものです。
現在の描画とは細部が異なる場合があります。

> 注意: 以下のスクリーンショットでは、インストールのセクションでオプションとしている「海岸線データ」「都市アウトラインデータ」「AKARI IR bands データ」を有効にしています。

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot1.png" alt="松江駅から見た夜空、Summer Triangle のアステリズム、画面右側の雲、都市アウトライン" width="100%" /></td>
    <td valign="top"><p><code>zstarview -p "Matsue Station" -A5 -Znw -a0.4 -P0.4</code> で表示した、日本の一都市の夜空です。恒星 Vega のあたりにマウスをホバーしているため、この星が属するアステリズムの Summer Triangle が表示されています。建物は <em>urban outline</em> として表示されています。</p><p>画面右側には雲が表示され、円形のドットを並べたハーフトーンで表現されています。雲は、GOES や Himawari などの衛星データから推定した雲量をもとに描画しています。</p><p>建物のデータは基本的には Overture Maps から取得しますが、日本の都市のいくつかについては、<a href="../README_ja-JP.md#plateau-building-data-preparation">PLATEAU のデータ</a>も利用可能です。このスクリーンショットでは PLATEAU のデータを利用しています。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot4.png" alt="アトランタ空港と上空の混雑した空域、航空機の軌跡と昇らない領域" width="100%" /></td>
    <td valign="top"><p>アトランタ空港とその上空の画面で、忙しい空域を飛ぶ10機以上の航空機が表示されています。航空機の軌跡は紫色のリボンで描かれ、地平線の下にある楕円は <em>never-rises</em>（昇らない領域）、つまり天球上で水平線より上に来ることがない領域を示しています。天の赤道と黄道は点線の基準線として表示されています。</p><p>天の赤道は、天の北極と天の南極を結ぶ地軸を法線とする大円です。never-rises の境界は、同じ地軸を中心軸とする小円です。</p><p>注意: 航空機をリボンとして表示するには、CLI では <code>-a 0.4</code> などを指定してください。<code>zstarview-gui</code> では、起動時のダイアログの <code>Overlays</code> タブで <code>Aircraft opacity</code> を 0.4 などに設定してください。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot3.png" alt="ユウニ塩湖から見た空とほぼ円形の水平線" width="100%" /></td>
    <td valign="top"><p>世界で最も平らな場所の一つとして知られているユウニ塩湖から、天頂を見上げた画面です。周囲をぐるっと囲む水平線は、地平線の高低差がほとんどないため、ほぼ円形に描かれています。<code>-V9</code> により、視等級 9 までの恒星を表示しています。また、<code>-s5</code> により、小さな恒星を目立たせるため、デフォルトより少し大きめに恒星を表示しています。</p><p>注意: 等級上限を大きくすると描画時間も増えます。<a href="cli-celestial-ja_JP.md#about-magnitude-limit">等級上限オプションについて</a>も参照してください。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot9.png" alt="地表から108mのタワーで高度-5度から見た夜間光の地面" width="100%" /></td>
    <td valign="top"><p>地表から 108m の高さにあるタワーから、視線方向の高度を少し下向き（<code>-A-5</code>、-5 度）にした画面です。Kobe Port Tower のように、いくつかのタワーは内部データベースに場所と高さが登録されています。利用可能なタワー名は <code>--list-viewpoints t</code> で、山名は <code>--list-viewpoints m</code> で一覧表示できます。この画面では、<em>night light</em>（宇宙から見た地表の明るさのデータ）によって地面が光る効果も付けています。これにより、都市の部分と海の部分で明るさが異なって見えます。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot11.png" alt="-A40で見上げた松江上空の左上に表示された月をマウスホバーで5倍に拡大" width="100%" /></td>
    <td valign="top"><p>視線方向の高度を変更するオプションに <code>-A40</code> を指定し、空を見上げた状況を示しています。FOV は画面の端で 90 度になるため、少し魚眼レンズのような効果が現れます。月は画面の左上側に表示されています。</p><p>月の上にマウスを移動すると、月が通常の見かけの大きさの 5 倍に拡大されます。通常の動作では、表示時刻に対応する月のシミュレーション画像を NASA Scientific Visualization Studio の Dial-A-Moon API から取得します。このシミュレーション画像には、影として月の位相が含まれています。</p><p>デフォルトの明るい天体モードでは、通常サイズの月は明るい側の外周弧と明暗境界の弧による月相アウトラインで表示されます。視等級に応じて描かれる恒星や惑星とは異なり、月は見かけの角直径を基準に円盤として描かれます。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot10.png" alt="オリーブ色の稜線とハーフトーンの雲オーバーレイ" width="100%" /></td>
    <td valign="top"><p>スイスの山岳地帯を表示した画面で、切り立った地形が水平線付近に見えています。太陽が地平線より上にあるため、空には色が付いています。稜線はオリーブ色で描かれています。</p><p>雲の描画については、ヨーロッパの大部分は GOES や Himawari の衛星のカバー範囲外なので、デフォルトでは雲は描かれません。<code>--geo-satellite true</code> を指定すると、静止衛星の雲画像を加工して表示する実験的なサポートを利用できます。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot6.png" alt="zstarview-export-image で生成した Torifune と大阪城の表示" width="100%" /></td>
    <td valign="top"><p>GUI アプリではなく <code>zstarview-export-image</code> を用いて出力した星空画像です。オブジェクトの検索オプション <code>--search "Torifune"</code> を利用して、小天体の位置を表示しています。右側には日本の建物、大阪城も見えています。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot12.png" alt="東京タワーから見た台風と予報降水量" width="100%" /></td>
    <td valign="top"><p>この画面では、東京近辺を移動する台風を表示しています。赤いマーカーが台風を、青い斜めの線がモデルによって予測された降水量を示します。雲などの観測に基づくレイヤーとは異なり、この予測は実際の降水と異なることがあります。予測降水量の多さは線の本数で表現しています。</p><p>降水量表示は、<code>--precipitation-opacity</code> に0以外の値を指定すると有効になります。初めて有効にするときは、Open-Meteo Free APIの利用条件への同意が必要です。</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot14.png" alt="ブダペストから見たGlobal Meteor Networkの流星軌跡" width="100%" /></td>
    <td valign="top"><p>この画面では、ブダペスト上空で観測された Global Meteor Network の流星軌跡を表示しています。軌跡は観測時に記録された方向に描画され、明るく先細りする流星の経路と、<code>-56h</code> や <code>-67h</code> のような経過時間ラベルが表示されます。ペルセウス座流星群の時期に撮影されたため、特に多くの軌跡が表示されています。</p><p>流星軌跡は、<code>--meteor-trails-opacity</code> に正の値を指定すると有効になります。<code>--meteor-trails-max-candidates</code> で表示する軌跡数を制限できます。</p></td>
  </tr>
</table>

<details>
<summary>その他のスクリーンショット</summary>

世界各地の都市アウトラインと地形地平線の例を示す画像です。これらの PNG は GUI のスクリーンショットではなく、`zstarview-export-image` の出力です。場所・時刻・視線方向のメタデータは PNG のテキストチャンクに埋め込まれており、`exiftool` などで確認できます。

<table>
  <tr>
    <td align="center" width="25%"><img src="images/screenshot-sydney.png" alt="サーキュラー・キー（シドニー）" width="100%" /></td>
    <td align="center" width="25%"><img src="images/screenshot-tokyotower.png" alt="東京タワー付近（東京）" width="100%" /></td>
    <td align="center" width="25%"><img src="images/screenshot-mtfuji.png" alt="山中湖から望む富士山" width="100%" /></td>
    <td align="center" width="25%"><img src="images/screenshot-izumo.png" alt="出雲大社（島根）" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>サーキュラー・キー（シドニー）</sub></td>
    <td align="center"><sub>東京タワー付近（東京）</sub></td>
    <td align="center"><sub>山中湖から望む富士山</sub></td>
    <td align="center"><sub>出雲大社（島根）</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/screenshot-taipei.png" alt="台北101（台湾）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-marinabay.png" alt="マリーナベイ（シンガポール）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-burjkhalifa-nc.png" alt="ブルジュ・ハリファ付近（ドバイ）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-manarola.png" alt="マナローラ（イタリア）" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>台北101（台湾）</sub></td>
    <td align="center"><sub>マリーナベイ（シンガポール）</sub></td>
    <td align="center"><sub>ブルジュ・ハリファ付近（ドバイ）</sub></td>
    <td align="center"><sub>マナローラ（イタリア）</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/screenshot-jungfraujoch.png" alt="ユングフラウヨッホ（Jungfraujoch、スイス）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-sagradafamilia.png" alt="サグラダファミリア（バルセロナ）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-westminsterbridge.png" alt="ウェストミンスター橋（ロンドン）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-uyuni.png" alt="ユウニ塩湖の別の画像（ボリビア）" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>ユングフラウヨッホ（Jungfraujoch、スイス）</sub></td>
    <td align="center"><sub>サグラダファミリア（バルセロナ）</sub></td>
    <td align="center"><sub>ウェストミンスター橋（ロンドン）</sub></td>
    <td align="center"><sub>ユウニ塩湖の別の画像（ボリビア）</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/screenshot-wizardisland.png" alt="ウィザード島（オレゴン）" width="100%" /></td>
    <td align="center"><img src="images/screenshot-bliss.png" alt="Bliss の丘（カリフォルニア）" width="100%" /></td>
    <td colspan="2"></td>
  </tr>
  <tr>
    <td align="center"><sub>ウィザード島（オレゴン）</sub></td>
    <td align="center"><sub>Bliss の丘（カリフォルニア）</sub></td>
    <td colspan="2"></td>
  </tr>
</table>

注記: 一部のスクリーンショットでは PLATEAU データを使用しています。詳細は
[PLATEAU 建物データの準備](../README_ja-JP.md#plateau-building-data-preparation) を参照してください。
</details>
