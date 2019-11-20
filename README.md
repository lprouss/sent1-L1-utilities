MATLAB utilities for Sentinel-1 Level 1 datasets
================================================

This repository contains several MATLAB functions to read and prepare/process Sentinel-1 Level 1 datasets.

These functions were developed primarily to work with single-look complex (SLC) images acquired in TOPS mode, i.e. Sentinel-1's Interferometric Wide Swath (IW) or Extra Wide Swath (EW) modes. Expected compatibility with the StripMap (SM) and Wave (WV) modes, and with ground-range detected (GRD) images, is indicated for each function.


## List of functions

__NOTE:__ These functions may not work with MATLAB versions older than R2016b.

+ `parseSafeManifest`: Extract useful information from the `manifest.safe` file of a Sentinel-1 dataset. This function should work with any Sentinel-1 Level 1 dataset.
+ `readSent1Data`: Read data in a Sentinel-1 measurement TIFF file. This function should work with any Sentinel-1 Level 1 image.
+ `applySent1Lut`: Apply a calibration LUT to a Sentinel-1 SLC image. This function should be compatible with SLC images in all modes. It has not been tested with GRD images.
+ `private/getTIFFinfo`: Extract useful parameters from a TIFF file, e.g. information used to read it, and perform basic checks. This function should work with any Sentinel-1 measurement TIFF file, but also for RADARSAT-2 images in the BigTIFF format.

Additional functions in the `dev` branch (under development or not fully tested):

+ `mergeSent1TopsBursts`: Stitch (merge) bursts in a TOPS subswath. Applicable only to IW and EW SLC images.
+ `topsDeramping`: Deramp TOPS data in azimuth. Applicable only to IW and EW SLC images.

Some of these functions use the [`xmlExtract` utility](https://github.com/lprouss/xmlExtract), which provides an (arguably) easier way to extract data from XML files, like Sentinel-1's annotation files. This function, as well as its associated sub-functions, should thus be present in MATLAB's search path.

