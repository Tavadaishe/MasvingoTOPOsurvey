/**********************************************************************
 * Masvingo – Hydrothermal Alteration Mapping with Harmonized
 * Landsat Sentinel-2 (HLS v2.0) in Google Earth Engine
 *
 * Products
 *   1. Iron Oxide Index (IOI)            = Red / Blue
 *   2. Clay (Hydroxyl) Alteration Index  = SWIR1 / SWIR2
 *   3. Quartz Index (QI)                 = TIRS1 / TIRS2   (HLSL30 only)
 *   4. Ferrous Iron Index (FeI, support) = SWIR2/NIR + Green/Red
 *   5. Hydrothermal Alteration Index (HAI) – weighted combination of the
 *      percentile-normalised indices, plus a 4-class alteration map.
 *
 * Paste into the GEE Code Editor (https://code.earthengine.google.com)
 * and press Run. Exports appear in the "Tasks" tab – click RUN on each.
 **********************************************************************/

// ============================ 1. SETTINGS ============================
// Drone survey footprint (from the Masvingo web map). The AOI is the
// survey rectangle buffered, because the survey itself is only ~1.5 km
// wide (~50 x 40 HLS pixels). Replace with your own asset if needed, e.g.
// var aoi = ee.FeatureCollection('users/<you>/masvingo_aoi').geometry();
var surveyArea = ee.Geometry.Rectangle([30.73700, -20.05722, 30.75073, -20.04561]);
var BUFFER_M   = 10000;                       // AOI buffer in metres
var aoi        = surveyArea.buffer(BUFFER_M).bounds();

var START_DATE   = '2022-01-01';
var END_DATE     = '2026-01-01';
var DRY_MONTHS   = [5, 10];                   // May–Oct: low vegetation, few clouds
var MAX_CLOUD    = 30;                        // scene-level cloud cover (%)
var NDVI_MAX     = 0.35;                      // mask denser vegetation (set to 1 to disable)
var EXPORT_SCALE = 30;                        // HLS native resolution (m)
var EXPORT_CRS   = 'EPSG:32736';              // WGS84 / UTM zone 36S
var EXPORT_FOLDER = 'Masvingo_HLS_Alteration';
var KELVIN       = 273.15;                    // HLSL30 TIRS bands are in °C

// HAI weights (IOI, Clay, Quartz). Quartz falls back to 0 where missing.
var W_IOI = 0.35, W_CLAY = 0.40, W_QI = 0.25;

Map.centerObject(aoi, 12);

// ===================== 2. HLS COLLECTIONS & MASKING ==================
// Fmask bits: 1 cloud, 2 adjacent cloud/shadow, 3 cloud shadow,
//             4 snow/ice, 5 water
function maskHLS(img) {
  var fmask = img.select('Fmask');
  var bad = fmask.bitwiseAnd(1 << 1)
    .or(fmask.bitwiseAnd(1 << 2))
    .or(fmask.bitwiseAnd(1 << 3))
    .or(fmask.bitwiseAnd(1 << 4))
    .or(fmask.bitwiseAnd(1 << 5));
  return img.updateMask(bad.eq(0));
}

// Harmonise band names between Landsat (L30) and Sentinel-2 (S30)
var COMMON = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

var l30 = ee.ImageCollection('NASA/HLS/HLSL30/v002')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.calendarRange(DRY_MONTHS[0], DRY_MONTHS[1], 'month'))
  .filter(ee.Filter.lt('CLOUD_COVERAGE', MAX_CLOUD))
  .map(maskHLS);

var s30 = ee.ImageCollection('NASA/HLS/HLSS30/v002')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.calendarRange(DRY_MONTHS[0], DRY_MONTHS[1], 'month'))
  .filter(ee.Filter.lt('CLOUD_COVERAGE', MAX_CLOUD))
  .map(maskHLS);

var l30Optical = l30.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'], COMMON);
var s30Optical = s30.select(['B2', 'B3', 'B4', 'B8A', 'B11', 'B12'], COMMON);

var hls = l30Optical.merge(s30Optical);
print('HLSL30 scenes:', l30.size(), 'HLSS30 scenes:', s30.size(),
      'Total HLS scenes:', hls.size());

// Median dry-season reflectance composite
var composite = hls.median().clip(aoi);

// Thermal composite (Landsat TIRS in HLSL30 only), converted to Kelvin
var thermal = l30.select(['B10', 'B11'], ['tir1', 'tir2'])
  .median().add(KELVIN).clip(aoi);

// ========================== 3. INDICES ===============================
var ndvi = composite.normalizedDifference(['nir', 'red']).rename('NDVI');
var vegMask = ndvi.lt(NDVI_MAX);

// Iron oxide index: hematite/goethite absorb strongly in the blue
var ioi = composite.select('red').divide(composite.select('blue')).rename('IOI');

// Clay / hydroxyl alteration index: Al-OH & Mg-OH absorb near 2.2 µm (SWIR2)
var clay = composite.select('swir1').divide(composite.select('swir2')).rename('CAI');

// Ferrous iron index (supporting indicator)
var ferrous = composite.select('swir2').divide(composite.select('nir'))
  .add(composite.select('green').divide(composite.select('red')))
  .rename('FeI');

// Quartz index from TIRS brightness temperature ratio. HLS has no
// ASTER-like 8–9.5 µm bands, so this is a proxy for silica-rich rocks.
var quartz = thermal.select('tir1').divide(thermal.select('tir2')).rename('QI');

var rawIndices = ee.Image.cat([ioi, clay, quartz, ferrous, ndvi])
  .updateMask(vegMask)
  .toFloat();

// ===================== 4. NORMALISATION & HAI ========================
// Robust 2–98 % percentile stretch to 0–1 so indices are comparable
function normalise(img, band) {
  var p = img.select(band).reduceRegion({
    reducer: ee.Reducer.percentile([2, 98]),
    geometry: aoi, scale: EXPORT_SCALE, maxPixels: 1e10, bestEffort: true
  });
  var lo = ee.Number(p.get(band + '_p2'));
  var hi = ee.Number(p.get(band + '_p98'));
  return img.select(band).subtract(lo).divide(hi.subtract(lo))
    .clamp(0, 1).rename(band + '_n');
}

var ioiN   = normalise(rawIndices, 'IOI');
var clayN  = normalise(rawIndices, 'CAI');
var qiN    = normalise(rawIndices, 'QI');
var ferN   = normalise(rawIndices, 'FeI');

// Where QI is unavailable, re-weight IOI and Clay only
var qiFilled = qiN.unmask(-1);
var hai = ee.Image().expression(
  'q < 0 ? (wi*i + wc*c) / (wi + wc) : wi*i + wc*c + wq*q', {
    i: ioiN, c: clayN, q: qiFilled, wi: W_IOI, wc: W_CLAY, wq: W_QI
  }).updateMask(vegMask).rename('HAI');

// Classify HAI by statistical thresholds (mean + kσ)
var haiStats = hai.reduceRegion({
  reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true),
  geometry: aoi, scale: EXPORT_SCALE, maxPixels: 1e10, bestEffort: true
});
var mu = ee.Number(haiStats.get('HAI_mean'));
var sd = ee.Number(haiStats.get('HAI_stdDev'));
print('HAI mean / std:', mu, sd);

var haiClass = ee.Image(1)
  .where(hai.gt(mu.add(sd)), 2)               // weak
  .where(hai.gt(mu.add(sd.multiply(1.5))), 3) // moderate
  .where(hai.gt(mu.add(sd.multiply(2))), 4)   // strong
  .updateMask(hai.mask())
  .rename('HAI_class')
  .toByte();

// Area per alteration class (ha)
var areaHa = ee.Image.pixelArea().divide(1e4).addBands(haiClass)
  .reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: aoi, scale: EXPORT_SCALE, maxPixels: 1e10, bestEffort: true
  });
print('Area (ha) per HAI class [1 background, 2 weak, 3 moderate, 4 strong]:', areaHa);

// ======================= 5. VISUALISATION ============================
var PAL_IRON   = ['#fff5eb', '#fdbe85', '#fd8d3c', '#d94701', '#7f2704'];
var PAL_CLAY   = ['#f7fcfd', '#ccece6', '#66c2a4', '#238b45', '#00441b'];
var PAL_QUARTZ = ['#fcfbfd', '#dadaeb', '#9e9ac8', '#6a51a3', '#3f007d'];
var PAL_HAI    = ['#313695', '#74add1', '#ffffbf', '#f46d43', '#a50026'];
var PAL_CLASS  = ['#d9d9d9', '#fee08b', '#fc8d59', '#d73027'];

Map.addLayer(composite, {bands: ['red', 'green', 'blue'], min: 0.02, max: 0.25},
  'HLS true colour', true);
Map.addLayer(composite, {bands: ['swir2', 'nir', 'red'], min: 0.02, max: 0.40},
  'HLS false colour (SWIR2-NIR-Red)', false);
// Sabins-style ratio composite: R = Clay, G = Iron oxide, B = Ferrous
Map.addLayer(ee.Image.cat([clayN, ioiN, ferN]), {min: 0, max: 1},
  'Ratio RGB (Clay / Iron oxide / Ferrous)', false);

Map.addLayer(ioiN,   {min: 0, max: 1, palette: PAL_IRON},   'Iron Oxide Index', false);
Map.addLayer(clayN,  {min: 0, max: 1, palette: PAL_CLAY},   'Clay Alteration Index', false);
Map.addLayer(qiN,    {min: 0, max: 1, palette: PAL_QUARTZ}, 'Quartz Index', false);
Map.addLayer(ferN,   {min: 0, max: 1, palette: PAL_IRON},   'Ferrous Iron Index', false);
Map.addLayer(hai,    {min: 0, max: 1, palette: PAL_HAI},    'Hydrothermal Alteration Index', true);
Map.addLayer(haiClass, {min: 1, max: 4, palette: PAL_CLASS}, 'Hydrothermal Alteration Classes', false);

Map.addLayer(ee.Image().paint(aoi, 0, 2), {palette: ['black']}, 'AOI');
Map.addLayer(ee.Image().paint(surveyArea, 0, 2), {palette: ['cyan']}, 'Drone survey area');

// ---- Legends ----
function colourBar(title, palette, lowLabel, highLabel) {
  var bar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '150x10', format: 'png',
             min: 0, max: 1, palette: palette},
    style: {stretch: 'horizontal', margin: '0 8px', maxHeight: '18px'}
  });
  var labels = ui.Panel([
    ui.Label(lowLabel, {margin: '2px 8px', fontSize: '11px'}),
    ui.Label(highLabel, {margin: '2px 8px', fontSize: '11px', textAlign: 'right',
                         stretch: 'horizontal'})
  ], ui.Panel.Layout.flow('horizontal'));
  return ui.Panel([ui.Label(title, {fontWeight: 'bold', fontSize: '12px',
                                    margin: '6px 8px 2px'}), bar, labels]);
}

function classLegend(title, names, colours) {
  var panel = ui.Panel([ui.Label(title, {fontWeight: 'bold', fontSize: '12px',
                                         margin: '6px 8px 2px'})]);
  names.forEach(function (name, i) {
    panel.add(ui.Panel([
      ui.Label('', {backgroundColor: colours[i], padding: '7px', margin: '2px 6px 2px 8px'}),
      ui.Label(name, {fontSize: '11px', margin: '3px 0'})
    ], ui.Panel.Layout.flow('horizontal')));
  });
  return panel;
}

var legend = ui.Panel({style: {position: 'bottom-left', padding: '6px', width: '210px'}});
legend.add(ui.Label('Masvingo – HLS Alteration', {fontWeight: 'bold', fontSize: '14px'}));
legend.add(colourBar('Iron Oxide Index (Red/Blue)', PAL_IRON, 'Low', 'High'));
legend.add(colourBar('Clay Index (SWIR1/SWIR2)', PAL_CLAY, 'Low', 'High'));
legend.add(colourBar('Quartz Index (TIR1/TIR2)', PAL_QUARTZ, 'Low', 'High'));
legend.add(colourBar('Hydrothermal Alteration Index', PAL_HAI, 'Low', 'High'));
legend.add(classLegend('Alteration Classes',
  ['Background', 'Weak (> μ+1σ)', 'Moderate (> μ+1.5σ)', 'Strong (> μ+2σ)'], PAL_CLASS));
Map.add(legend);

// ===================== 6. HISTOGRAMS (Console) =======================
print(ui.Chart.image.histogram({
  image: ee.Image.cat([ioiN, clayN, qiN, hai]),
  region: aoi, scale: EXPORT_SCALE, maxBuckets: 100, maxPixels: 1e9
}).setOptions({title: 'Normalised index distributions',
               hAxis: {title: 'Normalised value (0–1)'}, vAxis: {title: 'Pixels'}}));

// ========================== 7. EXPORTS ===============================
function exportToDrive(image, name) {
  Export.image.toDrive({
    image: image,
    description: name,
    folder: EXPORT_FOLDER,
    fileNamePrefix: name,
    region: aoi,
    scale: EXPORT_SCALE,
    crs: EXPORT_CRS,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true, noData: -9999}
  });
}

// Individual float rasters (raw ratio values)
exportToDrive(ioi.updateMask(vegMask).toFloat(),     'Masvingo_HLS_IronOxideIndex');
exportToDrive(clay.updateMask(vegMask).toFloat(),    'Masvingo_HLS_ClayAlterationIndex');
exportToDrive(quartz.updateMask(vegMask).toFloat(),  'Masvingo_HLS_QuartzIndex');
exportToDrive(ferrous.updateMask(vegMask).toFloat(), 'Masvingo_HLS_FerrousIronIndex');
exportToDrive(hai.toFloat(),                         'Masvingo_HLS_HydrothermalAlterationIndex');

// Alteration classes (Byte: 1 background, 2 weak, 3 moderate, 4 strong)
Export.image.toDrive({
  image: haiClass, description: 'Masvingo_HLS_HydrothermalAlterationClasses',
  folder: EXPORT_FOLDER, fileNamePrefix: 'Masvingo_HLS_HydrothermalAlterationClasses',
  region: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13,
  fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true, noData: 0}
});

// Multi-band stack: raw + normalised indices + HAI (+ class as float)
exportToDrive(ee.Image.cat([rawIndices, ioiN, clayN, qiN, ferN, hai,
                            haiClass.toFloat()]).toFloat(),
              'Masvingo_HLS_AlterationStack');

// Surface reflectance composite for reference
exportToDrive(composite.toFloat(), 'Masvingo_HLS_DrySeasonComposite');

// Styled RGB maps (ready for GIS / web map overlays)
exportToDrive(ioiN.visualize({min: 0, max: 1, palette: PAL_IRON}),     'Masvingo_HLS_IronOxide_RGB');
exportToDrive(clayN.visualize({min: 0, max: 1, palette: PAL_CLAY}),    'Masvingo_HLS_Clay_RGB');
exportToDrive(qiN.visualize({min: 0, max: 1, palette: PAL_QUARTZ}),    'Masvingo_HLS_Quartz_RGB');
exportToDrive(hai.visualize({min: 0, max: 1, palette: PAL_HAI}),       'Masvingo_HLS_HAI_RGB');
exportToDrive(haiClass.visualize({min: 1, max: 4, palette: PAL_CLASS}), 'Masvingo_HLS_HAIClasses_RGB');

// Optional: save to your GEE assets instead of Drive (uncomment)
// Export.image.toAsset({
//   image: ee.Image.cat([ioiN, clayN, qiN, hai, haiClass.toFloat()]),
//   description: 'Masvingo_HLS_Alteration_Asset',
//   assetId: 'users/<your_username>/Masvingo_HLS_Alteration',
//   region: aoi, scale: EXPORT_SCALE, crs: EXPORT_CRS, maxPixels: 1e13
// });

// AOI boundary as shapefile
Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(aoi, {name: 'Masvingo_AOI'})]),
  description: 'Masvingo_HLS_AOI', folder: EXPORT_FOLDER, fileFormat: 'SHP'
});
