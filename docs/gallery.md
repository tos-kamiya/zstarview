# zstarview Screenshot Gallery

These screenshots were taken with earlier versions of zstarview, so some
details of the current rendering may differ.

> Note: The screenshots below have the optional coastline data and
> urban outline data and AKARI IR bands data described in the Installation
> section enabled.

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot1.png" alt="Night sky over Matsue Station with the Summer Triangle asterism, clouds on the right, and urban outline" width="100%" /></td>
    <td valign="top"><p>This view shows the night sky over a Japanese city, displayed with <code>zstarview -p "Matsue Station" -A5 -Znw -a0.4 -P0.4</code>. The mouse is hovering near Vega, so the asterism it belongs to, the Summer Triangle, is shown. Buildings are displayed as an <em>urban outline</em>.</p><p>Clouds are visible on the right side of the image, rendered as a halftone pattern of circular dots. Clouds are rendered from cloud amounts estimated from satellite data such as GOES and Himawari.</p><p>Building data is generally obtained from Overture Maps, but <a href="../README.md#plateau-building-data-preparation">PLATEAU data</a> is also available for some Japanese cities. This screenshot uses PLATEAU data.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot4.png" alt="Atlanta airport and busy airspace with aircraft trails and the never-rises region" width="100%" /></td>
    <td valign="top"><p>This view shows Atlanta airport and the busy airspace above it, with more than ten aircraft visible. Aircraft trails are rendered as purple ribbons, and the ellipse below the horizon marks the <em>never-rises</em> region: the part of the celestial sphere that never comes above the horizon. The celestial equator and ecliptic are shown as dashed reference lines.</p><p>The celestial equator is a great circle whose normal axis is the Earth's rotation axis, the line connecting the north and south celestial poles. The never-rises boundary is a smaller circle around the same axis.</p><p>Note: To display aircraft as ribbons, specify <code>-a 0.4</code> or another opacity value in the CLI. In <code>zstarview-gui</code>, set <code>Aircraft opacity</code> to 0.4 or another value in the <code>Atmosphere</code> tab of the startup dialog.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot3.png" alt="View of the sky and nearly circular horizon from Salar de Uyuni" width="100%" /></td>
    <td valign="top"><p>This view looks straight up from Salar de Uyuni, known as one of the flattest places in the world. The horizon forms an almost complete circle around the view because the difference in elevation along the horizon is very small. The <code>-V9</code> option displays stars down to visual magnitude 9, while <code>-s5</code> makes the stars slightly larger than the default so that the smaller stars stand out.</p><p>Note: higher magnitude limits increase rendering time. See <a href="cli-celestial.md#about-magnitude-limit">About magnitude limit</a>.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot9.png" alt="Night-light ground glow viewed from a 108-meter tower at -5 degrees altitude" width="100%" /></td>
    <td valign="top"><p>This view is from a tower 108 meters above the ground, with the viewing altitude lowered slightly to <code>-A-5</code> degrees. Some towers, such as Kobe Port Tower, have their locations and heights registered in the internal database. Use <code>--list-viewpoints t</code> to list the available tower names, or <code>--list-viewpoints m</code> to list the available mountain names. This view also uses <em>night lights</em>, a dataset of the brightness of the ground as seen from space, which makes the ground glow differently over urban and sea areas.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot11.png" alt="Moon enlarged by mouse hover in the upper-left of an upward-looking view over Matsue" width="100%" /></td>
    <td valign="top"><p>This view uses <code>-A40</code> to change the altitude of the viewing direction and shows a view looking upward into the sky. The field of view reaches 90 degrees at the edge of the screen, producing a subtle fisheye-lens effect. The Moon is displayed in the upper-left part of the image.</p><p>The Moon is enlarged to 5x its normal apparent size after the mouse is moved over it. During normal operation, a simulated Moon image for the display time is fetched from NASA Scientific Visualization Studio's Dial-A-Moon API. The simulation includes the Moon's phase as the shadowed portion of the image.</p><p>In the default bright-body mode, the normal Moon marker uses a compact phase-aware outline. Unlike stars and planets, which are drawn according to visual magnitude, the Moon is rendered as a disk using its apparent angular size.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot10.png" alt="Mountain terrain, colored sky, and halftone cloud overlay in Switzerland" width="100%" /></td>
    <td valign="top"><p>This view shows a mountainous region of Switzerland, where steep terrain is visible along the horizon. Because the Sun is above the horizon, the sky is colored rather than fully dark. Terrain ridges are drawn in olive.</p><p>For cloud rendering, most of Europe is outside the coverage of GOES and Himawari, so clouds are not drawn by default. Experimental support is available through <code>--geo-satellite true</code>, which processes cloud imagery from a geostationary satellite.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot6.png" alt="zstarview-export-image output showing Torifune and Osaka Castle" width="100%" /></td>
    <td valign="top"><p>This view is a star-field image generated with <code>zstarview-export-image</code>, rather than a GUI application screenshot. The object search option <code>--search "Torifune"</code> is used to show the position of the minor body. Osaka Castle, a Japanese building, is visible on the right.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot12.png" alt="Typhoon and forecast precipitation viewed from Tokyo Tower" width="100%" /></td>
    <td valign="top"><p>This view shows a typhoon moving near Tokyo. The red marker shows the typhoon, and the blue diagonal lines show precipitation predicted by a forecast model. Unlike observation-based layers such as clouds, this forecast may differ from the precipitation that actually occurs. The number of lines represents the forecast amount of precipitation.</p><p>Forecast precipitation is enabled by setting <code>--precipitation-opacity</code> to a positive value. The first time it is enabled, you must agree to the Open-Meteo Free API terms.</p></td>
  </tr>
</table>

<table>
  <tr>
    <td valign="top" width="33%"><img src="images/screenshot14.png" alt="Global Meteor Network trails viewed from Budapest" width="100%" /></td>
    <td valign="top"><p>This view shows observed Global Meteor Network meteor trails over Budapest. The trails are drawn in the directions recorded at observation time, with brighter tapered marks for the meteor paths and compact age labels such as <code>-56h</code> and <code>-67h</code>. The image was captured during the Perseid meteor shower, so an especially large number of trails is visible.</p><p>Meteor trails are enabled with a positive <code>--meteor-trails-opacity</code> value. Use <code>--meteor-trails-max-candidates</code> to limit the number of displayed trails.</p></td>
  </tr>
</table>

<details>
<summary>Other screenshots</summary>

These screenshots show urban outline and terrain horizon examples from several places worldwide. They are `zstarview-export-image` outputs, not GUI screenshots, and they carry location/time/view metadata in embedded PNG text chunks that tools such as `exiftool` can inspect:

<table>
  <tr>
    <td align="center" width="25%"><img src="images/screenshot-sydney.png" alt="Circular Quay, Sydney" width="100%" /></td>
    <td align="center" width="25%"><img src="images/screenshot-tokyotower.png" alt="Near Tokyo Tower, Tokyo" width="100%" /></td>
    <td align="center" width="25%"><img src="images/screenshot-mtfuji.png" alt="View of Mt. Fuji, Japan" width="100%" /></td>
    <td align="center" width="25%"><img src="images/screenshot-izumo.png" alt="Izumo Taisha, Japan" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Circular Quay, Sydney</sub></td>
    <td align="center"><sub>Near Tokyo Tower, Tokyo</sub></td>
    <td align="center"><sub>View of Mt. Fuji, Japan</sub></td>
    <td align="center"><sub>Izumo Taisha, Japan</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/screenshot-taipei.png" alt="Taipei 101, Taiwan" width="100%" /></td>
    <td align="center"><img src="images/screenshot-marinabay.png" alt="Marina Bay, Singapore" width="100%" /></td>
    <td align="center"><img src="images/screenshot-burjkhalifa-nc.png" alt="Near Burj Khalifa, Dubai" width="100%" /></td>
    <td align="center"><img src="images/screenshot-manarola.png" alt="Manarola, Italy" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Taipei 101, Taiwan</sub></td>
    <td align="center"><sub>Marina Bay, Singapore</sub></td>
    <td align="center"><sub>Near Burj Khalifa, Dubai</sub></td>
    <td align="center"><sub>Manarola, Italy</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/screenshot-jungfraujoch.png" alt="Jungfraujoch, Switzerland" width="100%" /></td>
    <td align="center"><img src="images/screenshot-sagradafamilia.png" alt="Sagrada Familia, Barcelona" width="100%" /></td>
    <td align="center"><img src="images/screenshot-westminsterbridge.png" alt="Westminster Bridge, London" width="100%" /></td>
    <td align="center"><img src="images/screenshot-uyuni.png" alt="Alternate view of Salar de Uyuni, Bolivia" width="100%" /></td>
  </tr>
  <tr>
    <td align="center"><sub>Jungfraujoch, Switzerland</sub></td>
    <td align="center"><sub>Sagrada Familia, Barcelona</sub></td>
    <td align="center"><sub>Westminster Bridge, London</sub></td>
    <td align="center"><sub>Alternate view of Salar de Uyuni, Bolivia</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/screenshot-wizardisland.png" alt="Wizard Island, Oregon" width="100%" /></td>
    <td align="center"><img src="images/screenshot-bliss.png" alt="Bliss hill, California" width="100%" /></td>
    <td colspan="2"></td>
  </tr>
  <tr>
    <td align="center"><sub>Wizard Island, Oregon</sub></td>
    <td align="center"><sub>Bliss hill, California</sub></td>
    <td colspan="2"></td>
  </tr>
</table>

Note: Some of these screenshots use PLATEAU data. See [PLATEAU Building Data Preparation](../README.md#plateau-building-data-preparation) for details.
</details>
