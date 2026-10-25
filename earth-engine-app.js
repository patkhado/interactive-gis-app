// ============================================================
// GIS-Based Multi-Criteria Suitability Modelling for Agriculture
// Waterberg District Municipality, Limpopo, South Africa
// Platform: Google Earth Engine
// ============================================================


 //1) Study area (AOI): Waterberg//


var waterberg = aoi
  .filter(ee.Filter.eq('GAUL1_NAME', 'Limpopo'))//Province filter//
  .filter(ee.Filter.eq('GAUL2_NAME', 'Waterberg'));//District filter//
Map.centerObject(waterberg, 8);
Map.addLayer(waterberg, {color: 'red'}, 'Waterberg');

//2) Climate: Rainfall (CHIRPS Daily -> annual totals -> mean 2005–2023)//

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY');
var rain = ee.ImageCollection.fromImages(ee.List.sequence(2005, 2023).map(function (y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  return chirps.filterDate(start, start.advance(1, 'year')).sum();
})).mean().clip(waterberg).rename('rain');
Map.addLayer(rain, {min: 300, max: 800, palette: ['white', 'blue']}, 'Rainfall');

//3) Climate: Temperature (ERA5-Land Monthly -> mean 2005–2023 -> °C)//

var temp = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterDate('2005-01-01', '2023-12-31')
  .select('temperature_2m').mean()
  .subtract(273.15).clip(waterberg).rename('temp');
Map.addLayer(temp, {min: 16, max: 24, palette: ['blue', 'yellow', 'red']}, 'Temperature');

 //4) Topography: Elevation (SRTM) + Slope//
 
var dem = ee.Image('USGS/SRTMGL1_003').clip(waterberg);
var slope = ee.Terrain.slope(dem).rename('slope');
Map.addLayer(dem, {min: 700, max: 2000}, 'Elevation', false);
Map.addLayer(slope, {min: 0, max: 20, palette: ['white', 'brown']}, 'Slope');

 //5) Soil: Texture (USDA texture class) + pH (H2O at 4C)//
 
var texture = ee.Image('OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02')
  .select('b0').clip(waterberg).rename('texture');
var ph = ee.Image('OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02')
  .select('b0').divide(10).clip(waterberg).rename('ph');
Map.addLayer(texture, {min: 1, max: 12, palette: ['brown', 'yellow', 'green']}, 'Texture', false);
Map.addLayer(ph, {min: 5, max: 8, palette: ['red', 'yellow', 'blue']}, 'pH');

//6) Land Cover: ESA WorldCover (first image of collection)//

var lc = ee.ImageCollection('ESA/WorldCover/v200').first().clip(waterberg).rename('lc');
Map.addLayer(lc, {}, 'Land cover', false);

//7) Water Features: Water occurrence + distance to water//

var occ = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').clip(waterberg);
var water = occ.gt(25).selfMask().unmask(0);
var distWater = water.fastDistanceTransform(512).sqrt()
  .multiply(ee.Image.pixelArea().sqrt()).rename('distWater').clip(waterberg);
Map.addLayer(distWater, {min: 0, max: 20000, palette: ['blue', 'white']}, 'Distance to water');

//8) Reclassification helper function//
  
function reclass(img, breaks, scores) {
  var out = ee.Image(scores[0]);
  for (var i = 0; i < breaks.length; i++) { out = out.where(img.gte(breaks[i]), scores[i + 1]); }
  return out;
}

// 9) Apply reclassification to each variable (scoring for suitability//

var rainS = reclass(rain, [350, 450, 600, 750], [1, 2, 3, 4, 5]);
var tempS = reclass(temp, [16, 19, 22, 25], [1, 3, 5, 4, 2]);
var slopeS = reclass(slope, [2, 5, 8, 15], [5, 4, 3, 2, 1]);

var phS = reclass(ph, [4.5, 5.5, 6.0, 7.5, 8.5], [1, 3, 4, 5, 4, 2]);
var texS = texture.remap([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [3, 3, 3, 4, 4, 5, 5, 5, 4, 3, 2, 1], 0);
var lcS = lc.remap([10, 20, 30, 40, 60], [2, 3, 4, 5, 1], 0);
var distS = reclass(distWater, [2000, 5000, 10000, 20000], [5, 4, 3, 2, 1]);

 //10)Constraints / Allowed area (protected areas + exclusions)//
 
Map.centerObject(waterberg, 8);

var protectedAreas = ee.FeatureCollection('WCMC/WDPA/current/polygons').filterBounds(waterberg);
var protectedMask = ee.Image(0).byte().paint(protectedAreas, 1).clip(waterberg);
var constraint = lc.neq(50).and(lc.neq(80)).and(lc.neq(90))
  .and(slope.lt(30)).and(protectedMask.eq(0));
Map.addLayer(constraint.selfMask(), {palette: ['green']}, 'Allowed area');

//11) AHP helper function: compute weights from comparison matrix//

function ahp(M) {
  var n = M.length;
  var colSum = [];
  var w = [];
  var lambda = 0;
  for (var a = 0; a < n; a++) {
    colSum[a] = 0;
    for (var b = 0; b < n; b++) {
      colSum[a] = colSum[a] + M[b][a];
    }
  }
  for (var c = 0; c < n; c++) {
    w[c] = 0;
    for (var d = 0; d < n; d++) {
      w[c] = w[c] + M[c][d] / colSum[d];
    }
    w[c] = w[c] / n;
  }
  for (var e = 0; e < n; e++) {
    var aw = 0;
    for (var f = 0; f < n; f++) {
      aw = aw + M[e][f] * w[f];
    }
    lambda = lambda + aw / w[e];
  }
  var RI = {3: 0.58, 4: 0.90, 5: 1.12, 6: 1.24}[n];
  var CI = (lambda / n - n) / (n - 1);
  print('Consistency ratio (need < 0.10)', CI / RI);
  return w;
}

//12) Define AHP matrices and compute weights//

var Mrainfed = [
  [1,   2,   2,   3,   4,   5],
  [1/2, 1,   2,   2,   3,   4],
  [1/2, 1/2, 1,   2,   2,   3],
  [1/3, 1/2, 1/2, 1,   2,   2],
  [1/4, 1/3, 1/2, 1/2, 1,   2],
  [1/5, 1/4, 1/3, 1/2, 1/2, 1]
];
var wRainfed = ahp(Mrainfed);
print('Rainfed weights', wRainfed);

var Mirrig = [
  [1,   2,   2,   3,   4],
  [1/2, 1,   2,   2,   3],
  [1/2, 1/2, 1,   2,   2],
  [1/3, 1/2, 1/2, 1,   2],
  [1/4, 1/3, 1/2, 1/2, 1]
];
var wIrrig = ahp(Mirrig);
print('Irrigated weights', wIrrig);

//13) Build weighted suitability scores (continuous 1-5 scale)//

var rainfedScore = rainS.multiply(wRainfed[0]).add(texS.multiply(wRainfed[1])).add(slopeS.multiply(wRainfed[2]))
  .add(tempS.multiply(wRainfed[3])).add(phS.multiply(wRainfed[4])).add(lcS.multiply(wRainfed[5]))
  .updateMask(constraint).clip(waterberg).rename('rainfed');
var irrigScore = distS.multiply(wIrrig[0]).add(slopeS.multiply(wIrrig[1])).add(texS.multiply(wIrrig[2]))
  .add(phS.multiply(wIrrig[3])).add(tempS.multiply(wIrrig[4])).updateMask(constraint).clip(waterberg).rename('irrigated');
Map.addLayer(rainfedScore, {min: 1, max: 5, palette: ['red', 'yellow', 'green']}, 'Rainfed score');
Map.addLayer(irrigScore, {min: 1, max: 5, palette: ['red', 'yellow', 'green']}, 'Irrigated score');

function classify(score) {
  return ee.Image(1).where(score.gte(2), 2).where(score.gte(3), 3).where(score.gte(4), 4)
    .updateMask(score.mask()).clip(waterberg).rename('class');
}
var rainfedClass = classify(rainfedScore);
var irrigClass = classify(irrigScore);
var classVis = {min: 1, max: 4, palette: ['d7191c', 'fdae61', 'a6d96a', '1a9641']};
Map.addLayer(rainfedClass, classVis, 'Rainfed classes');
Map.addLayer(irrigClass, classVis, 'Irrigated classes');

// 14) Classify continuous scores into discrete classes (1..4)//

function areaByClass(img, label) {
  var stats = ee.Image.pixelArea().divide(10000).addBands(img).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: waterberg.geometry(), scale: 250, maxPixels: 1e10, bestEffort: true
  });
  print(label, stats.get('groups'));
}
areaByClass(rainfedClass, 'Rainfed area (ha) by class');
areaByClass(irrigClass, 'Irrigated area (ha) by class');

 //15) Area statistics by class//
 
Export.image.toDrive({image: rainfedClass.toByte(), description: 'Waterberg_Rainfed_Suitability',
  folder: 'GEE_Waterberg', region: waterberg.geometry(), scale: 100, crs: 'EPSG:32735', maxPixels: 1e10});
Export.image.toDrive({image: irrigClass.toByte(), description: 'Waterberg_Irrigated_Suitability',
  folder: 'GEE_Waterberg', region: waterberg.geometry(), scale: 100, crs: 'EPSG:32735', maxPixels: 1e10});
  
  //16) Cross-check: cropland within suitability classes//
  
  function cropCheck(img, label) {
  var s = ee.Image.pixelArea().divide(10000).updateMask(lc.eq(40)).addBands(img).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: waterberg.geometry(), scale: 250, maxPixels: 1e10, bestEffort: true});
  print(label, s.get('groups'));
}
cropCheck(rainfedClass, 'Existing cropland (ha) by rainfed class');
cropCheck(irrigClass, 'Existing cropland (ha) by irrigated class');

 // 17) Sensitivity baseline: equal weights comparison//
 
var eqScore = rainS.add(texS).add(slopeS).add(tempS).add(phS).add(lcS).divide(6)
  .updateMask(constraint).clip(waterberg);
var eqClass = classify(eqScore);
Map.addLayer(eqClass, classVis, 'Rainfed classes (equal weights)', false);
areaByClass(eqClass, 'Rainfed area (ha), equal weights');

//18) Sensitivity analysis: nudge one weight and recompute classes//

function nudge(f) {
  var w = wRainfed.slice();
  w[0] = w[0] * f;
  var total = w[0] + w[1] + w[2] + w[3] + w[4] + w[5];
  return w.map(function (x) { return x / total; });
}
function rainfedClassFor(w) {
  var s = rainS.multiply(w[0]).add(texS.multiply(w[1])).add(slopeS.multiply(w[2]))
    .add(tempS.multiply(w[3])).add(phS.multiply(w[4])).add(lcS.multiply(w[5]))
    .updateMask(constraint).clip(waterberg);
  return classify(s);
}
areaByClass(rainfedClassFor(nudge(0.8)), 'Rainfall weight -20%');
areaByClass(rainfedClassFor(nudge(1.2)), 'Rainfall weight +20%');

function nudgeI(f) {
  var w = wIrrig.slice();
  w[0] = w[0] * f;
  var total = w[0] + w[1] + w[2] + w[3] + w[4];
  return w.map(function (x) { return x / total; });
}
function irrigClassFor(w) {
  var s = distS.multiply(w[0]).add(slopeS.multiply(w[1])).add(texS.multiply(w[2]))
    .add(phS.multiply(w[3])).add(tempS.multiply(w[4]))
    .updateMask(constraint).clip(waterberg);
  return classify(s);
}
areaByClass(irrigClassFor(nudgeI(0.8)), 'Irrigated: water distance weight -20%');
areaByClass(irrigClassFor(nudgeI(1.2)), 'Irrigated: water distance weight +20%');

//19) Export extra rasters (scores & factor suitability layers)//
function exp(img, name) {
  Export.image.toDrive({image: img.toFloat().clip(waterberg), description: name, folder: 'GEE_Waterberg',
    region: waterberg.geometry(), scale: 100, crs: 'EPSG:32735', maxPixels: 1e10});
}
exp(rainfedScore, 'Rainfed_Score'); exp(irrigScore, 'Irrigated_Score');
exp(rainS, 'S_Rainfall'); exp(tempS, 'S_Temperature'); exp(slopeS, 'S_Slope');
exp(phS, 'S_pH'); exp(texS, 'S_Texture'); exp(lcS, 'S_LandCover'); exp(distS, 'S_DistWater');
Export.table.toDrive({collection: waterberg, description: 'Waterberg_Boundary', folder: 'GEE_Waterberg', fileFormat: 'SHP'});
