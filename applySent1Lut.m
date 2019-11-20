function imgCal = applySent1Lut( img, calAnnFile, lut, varargin )
% Apply a look-up table (LUT) to radiometrically calibrate a Sentinel-1 Level 1
% image.
%
% Inputs:
%   - img: Sentinel-1 Level 1 image with size [Nrg Naz], where 'Nrg' is the
%       number of range samples (pixels) and 'Naz' is the number of azimuth
%       lines.
%   - calAnnFile: calibration annotation file (XML) for the input Sentinel-1
%       image, with full path.
%   - lut: name of the desired LUT (case-sensitive), see Additional information.
%   - imgLut (optional): name of the LUT already applied to the image
%       (case-sensitive), see Additional information. If not provided, it is
%       assumed that no LUT has been applied to the image, i.e. it has the
%       original application-specific LUT applied during signal processing.
%       In this case, 'imgLut' is set to 'original'.
%   - verbose (optional): logical flag to display progress information, false
%       by default.
%
% Outputs:
%   - imgCal: calibrated Sentinel-1 with the desired LUT.
%
% Required functions (toolboxes and/or user-defined):
%   - xmlExtract (see https://github.com/lprouss/xmlExtract)
%
% Additional information:
%   This function uses the azimuth line and range sample (pixel) indexes in the
%   calibration annotation file to generate the LUT for the image. Therefore,
%   it only works if the complete image found in a measurement data file is
%   provided as input. The size of the input image is not validated.
%
%   For Sentinel-1, four calibration LUTs are available:
%   - 'dn': digital number for image pixels, without the application-specific
%       LUT applied during SAR processing;
%   - 'betaNought': backscattering coefficient in slant range geometry, also
%       known as radar brightness;
%   - 'sigmaNought': backscattering coefficient in ground range geometry;
%   - 'gamma': backscattering coefficient in the plane perpendicular to the
%       line-of-sight of the SAR system (also called gamma nought in the
%       literature).
%   When 'imgLut' is provided, this LUT is changed for the desired one.
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: November 2017
% Updated: January 2018, November 2019 (new version of the function based on
%   image pixel indexes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: add support for image subsets (regions of interest) by optionally
%   providing line/pixel offsets in azimuth and range as input?

%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% initialization

% validate number of inputs
narginchk( 3, 5 );

% default values for optional inputs
imgLut = 'original'; % image has original LUT
verbose = false; % no verbose

% size of input image array
validateattributes( img, {'numeric'}, {'2d'}, '', 'img', 1 );
[Nrg, Naz] = size( img );

% make sure the calibration annotation file exists
assert( exist( calAnnFile, 'file' ) == 2, ['Calibration annotation file ' ...
    'provided as input does not exist: %s.'], calAnnFile );

% make sure the desired LUT is valid
validLUT = {'sigmaNought', 'betaNought', 'gamma', 'dn'};
lut = validatestring( lut, validLUT, '', 'lut', 3 );

% process provided optional inputs, if any
for ni = 1:length( varargin )
    % assign input to correct variable using its data type
    if ischar( varargin{ni} )
        % string: image LUT
        validLUT = [validLUT, {'original'}];
        imgLut = validatestring( varargin{ni}, validLUT, '', 'imgLut', ni + 3 );
    elseif islogical( varargin{ni} )
        % logical: verbose flag
        verbose = varargin{ni};
    else
        % invalid type: return error
        error( ['Optional input #%d has an invalid type. Please refer to the ' ...
            'help text of this function for a description of inputs.'], ni );
    end
end

% return to caller is the requested LUT is already applied to the image
if strcmpi( lut, imgLut )
    fprintf( ['The requested LUT is already applied to the input image: ' ...
        'nothing to do.\n'] );
    return;
end


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% parse calibration annotation file and validate extracted information

% display progress information, if requested
if verbose
    fprintf( 'Extracting calibration data from the annotation file... ' );
end

% parsing information structure for the calibration annotation file
pinfo.tag = {'calibration', 'calibrationVector', 'line', 'pixel', lut};
pinfo.type = {'root', 'list', 'dblArr', 'dblArr', 'dblArr'};
pinfo.level = [0, 1, 2, 2, 2];
if ~strcmp( imgLut, 'original' )
    % LUT applied to input image: extract it too
    pinfo.tag = cat( 2, pinfo.tag, imgLut );
    pinfo.type = cat( 2, pinfo.type, 'dblArr' );
    pinfo.level = cat( 2, pinfo.level, 2 );
end
pinfo.dateStrFmt = 'yyyy-mm-ddTHH:MM:SS.FFFFFF';

% extract calibration data from the annotation file
pdata = xmlExtract( calAnnFile, pinfo );
Ncal = length( pdata.calibrationVector ); % number of cal. vectors

% get first vector of range pixels from calibration data
calRg = pdata.calibrationVector(1).pixel; % first range pixel vector
NrgCal = length( calRg ); % number of calibration range pixels

% make sure all vectors of calibration range pixels are identical
calRgAll = [pdata.calibrationVector.pixel]; % all range pixel vectors
assert( length( calRgAll ) / NrgCal == Ncal, ['The number of calibration ' ...
    'range pixels is not identical for all calibration vectors. This ' ...
    'function cannot proceed.'] );
calRgAll = unique( reshape( calRgAll, [NrgCal, Ncal] ).', 'rows' );
assert( size( calRgAll, 1 ) == 1, ['Calibration range pixels are not ' ...
    'identical for all calibration vectors. This function cannot proceed.'] );
clear calRgAll;

% vector of azimuth lines in calibration data
calAz = [pdata.calibrationVector.line];

% validate the number of range pixels and azimuth lines in calibration data
assert( calRg(end) + 1 == Nrg, ['Vector of calibration range pixels is ' ...
    'inconsistent with the input image: maximum pixel index is %d while ' ...
    'there are %d range pixels in the image.'], calRg(end) + 1, Nrg );
assert( calAz(end) + 1 >= Naz, ['Vector of calibration azimuth lines is ' ...
    'inconsistent with the input image: maximum line index is %d while ' ...
    'there are %d azimuth lines in the image.'], calAz(end) + 1, Naz );

% construct the desired LUT matrix
lutArray = [pdata.calibrationVector.(lut)];
assert( length( lutArray ) / NrgCal == Ncal, ['The number of calibration ' ...
    'values for the desired LUT is not identical for all calibration ' ...
    'vectors. This function cannot proceed.'] );
lutArray = reshape( lutArray, NrgCal, [] );

% construct the image LUT matrix, if necesary
if ~strcmp( imgLut, 'original' )
    imgLutArray = [pdata.calibrationVector.(imgLut)];
    assert( length( imgLutArray ) / NrgCal == Ncal, ['The number of ' ...
        'calibration values for the image LUT is not identical for all ' ...
        'calibration vectors. This function cannot proceed.'] );
    imgLutArray = reshape( imgLutArray, NrgCal, [] );
end
if verbose
    fprintf( 'Done!\n' );
end


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% interpolate LUT data

% display progress information, if requested
if verbose
    fprintf( 'Interpolating LUT array(s) to image... ' );
end

% create a mesh grid of calibration range pixels and azimuth lines
[calAzGrid, calRgGrid] = meshgrid( calAz, calRg );
clear calAz calRg;

% create a mesh grid of interpolation range pixels and azimuth lines
[imgAz, imgRg] = meshgrid( ( 0:Naz-1 ), ( 0:Nrg-1 ) );

% perform a bilinear interpolation of the desired LUT matrix
lutArrayInt = interp2( calAzGrid, calRgGrid, lutArray, imgAz, imgRg, ...
    'linear', NaN );
lutArrayInt = reshape( lutArrayInt, Nrg, Naz );
clear lutArray;

% perform a bilinear interpolation of the image LUT matrix, if necessary
if ~strcmp( imgLut, 'original' )
    imgLutArrayInt = interp2( calAzGrid, calRgGrid, imgLutArray, imgAz, ...
        imgRg, 'linear', NaN );
    imgLutArrayInt = reshape( imgLutArrayInt, Nrg, Naz );
    clear imgLutArray;
end
clear calAzGrid calRgGrid imgAz imgRg;
if verbose
    fprintf( 'Done!\n' );
end


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% apply LUT to input image

% display progress information, if requested
if verbose
    fprintf( 'Applying interpolated LUT to input image... ' );
end

% if necessary, extract phase information from input image array
if ~isreal( img )
    imgPhase = angle( img );
end

% calculate the magnitude squared of the image pixel values
img = abs( img ).^2;

% remove image LUT, if necessary
if ~strcmp( imgLut, 'original' )
    img = img .* imgLutArrayInt.^2;
    clear imgLutArrayInt;
end

% apply calibration LUT to image
imgCal = sqrt( img ./ lutArrayInt.^2 );
clear lutArrayInt img;

% reapply phase information to calibration image data
if exist( 'imgPhase', 'var' )
    imgCal = imgCal .* exp( 1j * imgPhase );
    clear imgPhase;
end
if verbose
    fprintf( 'Done!\n' );
end

