/**
 * Sentinel-2 MSI – band selection and export for mineral surveying
 * Google Earth Engine Code Editor script (https://code.earthengine.google.com)
 *
 * Area: Masvingo topographic survey site, Zimbabwe (buffered for regional context).
 *
 * Bands kept (Sentinel-2 L2A surface reflectance):
 *   B2  Blue        490 nm  10 m  – ferric iron absorption (with B4)
 *   B3  Green       560 nm  10 m  – ferrous iron index
 *   B4  Red         665 nm  10 m  – ferric iron / iron oxides (hematite, goethite)
 *   B5  Red Edge 1  705 nm  20 m  – iron features, vegetation separation
 *   B6  Red Edge 2  740 nm  20 m
 *   B7  Red Edge 3  783 nm  20 m
 *   B8  NIR         842 nm  10 m  – ferrous iron (Fe2+) broad absorption
 *   B8A NIR narrow  865 nm  20 m  – iron oxide crystal-field absorption
 *   B11 SWIR 1     1610 nm  20 m  – high reflectance for clays / carbonates
 *   B12 SWIR 2     2190 nm  20 m  – Al-OH, Mg-OH, CO3 absorption (clays, micas, carbonates)
 * Excluded: B1 (aerosol, 60 m), B9 (water vapour, 60 m), B10 (cirrus, not in L2A).
 */

// ---------------------------------------------------------------------------
// 1. Parameters
// ---------------------------------------------------------------------------
// Survey extent from the web map (index.html SURVEY_BOUNDS), buffered by 5 km.
// Replace `aoi` with your own drawn geometry/asset if needed.
var surveyExtent = ee.Geometry.Rectangle([30.73700, -20.05722, 30.75073, -20.04561]);
var aoi = surveyExtent.buffer(5000).bounds();

// Dry season (May–Oct) = minimal green vegetation and clouds in Zimbabwe,
// so bare rock/soil spectra are best exposed.
var START_DATE = '2025-05-01';
var END_DATE   = '2025-10-31';
var CLEAR_THRESHOLD = 0.60;   // Cloud Score+ cs_cdf threshold (higher = stricter)

var MINERAL_BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];

var EXPORT_FOLDER = 'GEE_Masvingo_Mineral';
var EXPORT_SCALE  = 10;        // metres (20 m bands are resampled to 10 m)
var EXPORT_CRS    = 'EPSG:32736'; // WGS 84 / UTM zone 36S (Masvingo)

// ---------------------------------------------------------------------------
// 2. Load Sentinel-2 SR + Cloud Score+ and mask clouds/shadows
// ---------------------------------------------------------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30));

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

function maskAndScale(img) {
  var clear = img.select('cs_cdf').gte(CLEAR_THRESHOLD);
  // Also drop water (6) and snow (11) using the Scene Classification Layer.
  var scl = img.select('SCL');
  var notWaterSnow = scl.neq(6).and(scl.neq(11));
  return ee.Image(img.select(MINERAL_BANDS)
    .divide(10000)                       // DN -> surface reflectance (0–1)
    .updateMask(clear.and(notWaterSnow))
    .copyProperties(img, ['system:time_start']));
}

var collection = s2
  .linkCollection(csPlus, ['cs_cdf'])
  .map(maskAndScale);

print('Scenes used:', collection.size());

// Median composite of the selected mineral bands, clipped to AOI
var composite = collection.median().select(MINERAL_BANDS).clip(aoi).toFloat();

// ---------------------------------------------------------------------------
// 3. Mineral / alteration indices (band ratios)
// ---------------------------------------------------------------------------
var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');

var indices = ee.Image.cat([
  // Ferric iron (Fe3+) – hematite, goethite, jarosite (gossans)
  composite.expression('B4 / B2', {B4: composite.select('B4'), B2: composite.select('B2')})
    .rename('FerricIron'),
  // Ferrous iron (Fe2+) – e.g. mafic/ultramafic rocks
  composite.expression('(B12 / B8) + (B3 / B4)', {
    B12: composite.select('B12'), B8: composite.select('B8'),
    B3: composite.select('B3'),   B4: composite.select('B4')
  }).rename('FerrousIron'),
  // Iron oxide ratio
  composite.expression('B11 / B8', {B11: composite.select('B11'), B8: composite.select('B8')})
    .rename('FerricOxides'),
  // Clay / hydroxyl (Al-OH, Mg-OH) and carbonate alteration – sericite, kaolinite
  composite.expression('B11 / B12', {B11: composite.select('B11'), B12: composite.select('B12')})
    .rename('ClayHydroxyl'),
  // Gossan index
  composite.expression('B11 / B4', {B11: composite.select('B11'), B4: composite.select('B4')})
    .rename('Gossan'),
  ndvi
]).toFloat();

// Mask dense vegetation (NDVI > 0.4) so the ratios reflect rock/soil
var bareMask = ndvi.lt(0.4);
var indicesBare = indices.updateMask(bareMask);

// ---------------------------------------------------------------------------
// 4. Visualisation
// ---------------------------------------------------------------------------
Map.centerObject(aoi, 13);
Map.addLayer(composite, {bands: ['B4', 'B3', 'B2'], min: 0.02, max: 0.30}, 'True colour (4-3-2)');
Map.addLayer(composite, {bands: ['B12', 'B11', 'B8A'], min: 0.05, max: 0.45}, 'SWIR geology (12-11-8A)', false);
Map.addLayer(composite, {bands: ['B12', 'B8', 'B2'], min: 0.02, max: 0.40}, 'Lithology (12-8-2)', false);
// Classic Sabins-style ratio RGB: R = clay, G = ferrous, B = ferric
Map.addLayer(indicesBare.select(['ClayHydroxyl', 'FerrousIron', 'FerricIron']),
  {min: [1.0, 1.5, 1.0], max: [1.6, 3.0, 2.5]}, 'Ratio RGB (Clay-Ferrous-Ferric)', false);
Map.addLayer(indicesBare.select('ClayHydroxyl'), {min: 1.0, max: 1.6, palette: ['blue', 'white', 'red']}, 'Clay / hydroxyl', false);
Map.addLayer(indicesBare.select('FerricIron'),   {min: 1.0, max: 2.5, palette: ['blue', 'white', 'red']}, 'Ferric iron', false);
Map.addLayer(ee.Image().paint(surveyExtent, 0, 2), {palette: ['yellow']}, 'Drone survey extent');
Map.addLayer(ee.Image().paint(aoi, 0, 2), {palette: ['red']}, 'AOI');

// ---------------------------------------------------------------------------
// 5. Exports (Tasks tab -> RUN)
// ---------------------------------------------------------------------------
// 5a. Selected Sentinel-2 mineral bands (10-band reflectance stack)
Export.image.toDrive({
  image: composite,
  description: 'S2_Mineral_Bands_Masvingo',
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'S2_mineral_bands_' + START_DATE + '_' + END_DATE,
  region: aoi,
  scale: EXPORT_SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});

// 5b. Mineral / alteration index stack (vegetation masked)
Export.image.toDrive({
  image: indicesBare,
  description: 'S2_Mineral_Indices_Masvingo',
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'S2_mineral_indices_' + START_DATE + '_' + END_DATE,
  region: aoi,
  scale: EXPORT_SCALE,
  crs: EXPORT_CRS,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF',
  formatOptions: {cloudOptimized: true}
});

// 5c. (Optional) Save the band stack as an Earth Engine asset for reuse.
// Replace the assetId with your own project path, then uncomment.
// Export.image.toAsset({
//   image: composite,
//   description: 'S2_Mineral_Bands_Masvingo_Asset',
//   assetId: 'projects/YOUR_PROJECT/assets/S2_mineral_bands_masvingo',
//   region: aoi,
//   scale: EXPORT_SCALE,
//   crs: EXPORT_CRS,
//   maxPixels: 1e13
// });
