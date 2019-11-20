Suggested workflow using the Sentinel-1 L1 utilities
====================================================

1) Extract useful information from SAFE manifest file for the dataset located in directory `S1dset`:

    `maniHead = parseSafeManifest( S1dset );`

2) Read image for the desired swath (`sidx`), polarization (`pidx`) and WV mode vignette index (`vidx`):

    `img = readSent1Data( fullfile( maniHead.rootDir, maniHead.dataFile{sidx,pidx,vidx} ), true );`

3) Apply calibration look-up table (LUT) `lut` to image:

    `imgCal = applySent1Lut( img, fullfile( maniHead.rootDir, maniHead.calAnn{sidx,pidx,vidx} ), lut, true );`

    If you decide to change the LUT to `lutNew`:

    `imgCalNew = applySent1Lut( imgCal, fullfile( maniHead.rootDir, maniHead.calAnn{sidx,pidx,vidx} ), lutNew, lut, true );`


