MATLAB utilities for Sentinel-1 Level 1 datasets
================================================

This repository contains several MATLAB functions to read and prepare/process Sentinel-1 Level 1 datasets. These functions were developed primarily to handle single-look complex (SLC) images acquired in TOPS mode, i.e. using Sentinel-1's Interferometric Wide Swath (IW) or Extra Wide Swath (EW) modes. Some functions may also work for the StripMap (SM) mode and for ground-range detected (GRD) images. However, support for the Wave (WV) mode is likely limited (never tested).

## List of functions

__NOTE:__ These functions may not work with MATLAB versions older than R2016b.

+ `parseSafeManifest`: Extract useful information from the `manifest.safe` file of a Sentinel-1 dataset. This function should work with any Sentinel-1 Level 1 dataset.
+ `readSent1Data`: Read data in a Sentinel-1 measurement TIFF file. As it is fairly general, this function should work with any Sentinel-1 Level image.
+ `applySent1Lut`: Apply a calibration LUT to a Sentinel-1 SLC image. This function works with TOPS images before and after bursts merging, and should also work with SM images.
+ `private/getTIFFinfo`: Extract useful parameters from a TIFF file, e.g. information used to read it, and perform basic checks. This function should work with any Sentinel-1 measurement TIFF file, but also for RADARSAT-2 images in the BigTIFF format.

Some of these functions use the [`xmlExtract` utility](https://github.com/lprouss/xmlExtract), which provides an (arguably) easier way to extract data from XML files, like Sentinel-1's annotation files. This function, as well as its associated sub-functions, should thus be present in MATLAB's search path.

