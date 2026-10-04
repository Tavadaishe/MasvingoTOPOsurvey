# MasvingoTOPOsurvey
Drone topographic survey web map – Masvingo

## Masvingo Project – Topological survey

Live web map: https://tavadaishe.github.io/MasvingoTOPOsurvey/

- `index.html` – Leaflet web map (orthomosaic always on, DSM/DTM and Gold Prospectivity overlays, Esri base maps)
- `tiles/ortho`, `tiles/dsm`, `tiles/dtm` – 512px XYZ tiles in WebP format
- `tiles/gold` – Gold Prospectivity Classes, 512px XYZ tiles in PNG format
- `gee/hls_alteration_indices.js` – Google Earth Engine script (Harmonized Landsat Sentinel-2, HLS v2.0) mapping the Iron Oxide, Clay Alteration, Quartz and Hydrothermal Alteration indices around the survey area, with GeoTIFF exports to Google Drive
