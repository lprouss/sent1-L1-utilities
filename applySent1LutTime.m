function imgCal = applySent1Lut( img, azTime, calAnnFile, lut, varargin )
% Apply a look-up table (LUT) to radiometrically calibrate a Sentinel-1 image.
%
% Inputs:
%   - img: Sentinel-1 image for one subswath (TOPS mode), with dimensions
%       [Nrg Naz] or [Nrg Naz Nb], where 'Nrg' and 'Naz' are respectively the
%       number of range and azimuth samples, and 'Nb' is the number of bursts
%       (TOPS image before bursts merging).
%   - azTime: structure containing azimuth timing information for the image,
%       see Additional information.
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
% Required functions (not part of MATLAB):
%   - xmlExtract (see https://github.com/lprouss/xmlExtract)
%
% Additional information:
%   The 'azTime' structure should contain the following parameters:
%   - 'T0': six-elements time vector for the first azimuth line in the image;
%   - 'trel': vector of azimuth time relative to 'T0', for each line in the
%       image. For a TOPS image before bursts merging, there should be one
%       vector per burst and 'trel' should have dimensions [Nb, Naz] ('Naz'
%       is the number of lines in each burst).
%
%   For Sentinel-1, four calibration LUTs are available:
%   - 'dn': (original) digital number for image pixels, without the
%       application-specific LUT applied during SAR processing;
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
% Updated: January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: support recovering the application-specific LUT, i.e. remove a
%   calibration LUT?
% TODO: support a header structure as input instead of an image (apply LUT to
%   each image in header)?

%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% initialization

% validate number of inputs
narginchk( 4, 6 );

% dimensions of input image array
dataDim = size( img );
Nrg = dataDim(1);
Naz = dataDim(2);

% number of bursts in the image
if length( dataDim ) == 3
    % TOPS image before bursts merging
    Nb = dataDim(3);
else
    % Stripmap SAR image or merged TOPS image
    Nb = 1;
end
clear dataDim;

% make sure the azimuth timing structure contains the required parameters
validTst = isstruct( azTime ) && all( isfield( azTime, {'T0', 'trel'} ) );
assert( validTst, ['The second input, azimuth timing information, is either ' ...
    'not a structure or does not contain the required parameters ("T0" and ' ...
    '"trel").'] );
T0 = azTime.T0(:).';
trel = azTime.trel;
clear azTime;

% validate the provided azimuth time information
validateattributes( T0, {'numeric'}, {'numel', 6, 'nonnegative'}, '', 'T0' );
validateattributes( trel, {'numeric'}, {'size', [Nb, Naz]}, '', 'trel' );

% make sure the calibration annotation file exists
assert( exist( calAnnFile, 'file' ) > 0, ['Calibration annotation file ' ...
    'provided as input does not exist: %s.'], calAnnFile );

% make sure the desired LUT is valid
validLUT = {'sigmaNought', 'betaNought', 'gamma', 'dn'};
lut = validatestring( lut, validLUT, '', 'lut', 3 );
clear validLUT;

% assign variable to the provided optional input(s)
imgLut = 'original'; % input image has original LUT
verbose = false; % no verbose
for ni = 1:nargin-4
    % assign current optional input to correct variable using its type
    if ischar( varargin{ni} )
        % string: LUT of the input image
        imgLut = varargin{ni};

        % make sure the image LUT is valid
        validLUT = {'sigmaNought', 'betaNought', 'gamma', 'dn', 'original'};
        imgLut = validatestring( imgLut, validLUT, '', 'imgLut', ni+4 );
        clear validLUT;
    elseif islogical( varargin{ni} )
        % logical: verbose flag
        verbose = varargin{ni};
    else
        % invalid type: return error
        error( ['Optional input #%d has an invalid type. Please refer to the ' ...
            'help text of this function for a description of inputs.'], ni );
    end
end


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% parse calibration annotation file and validate extracted information

% display progress information, if requested
if verbose
    fprintf( ['Extracting calibration information from the annotation file and ' ...
        'interpolating to image''s dimensions... '] );
end

% parsing information structure for the calibration annotation file
pinfo.tag = {'calibration', 'calibrationVector', 'azimuthTime', 'pixel', lut};
pinfo.type = {'root', 'list', 'dateStr', 'dblArr', 'dblArr'};
pinfo.level = [0, 1, 2, 2, 2];
if ~strcmp( imgLut, 'original' )
    % LUT applied to input image: extract it too
    pinfo.tag = cat( 2, pinfo.tag, imgLut );
    pinfo.type = cat( 2, pinfo.type, 'dblArr' );
    pinfo.level = cat( 2, pinfo.level, 2 );
end
pinfo.dateStrFmt = 'yyyy-mm-ddTHH:MM:SS.FFFFFF';

% extract information from calibration annotation file
pdata = xmlExtract( calAnnFile, pinfo );
Ncal = length( pdata.calibrationVector ); % number of cal. vectors
clear pinfo;

% get first vector of calibration range pixels from information structure
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

% make sure extract calibration data matches the image array
assert( calRg(end) + 1 == Nrg, ['Vector of calibration range pixels is ' ...
    'inconsistent with the input image array: maximum pixel value is %d ' ...
    'while there are %d range pixels in the data array.'], calRg(end) + 1, ...
    Nrg );

% construct the desired LUT matrix
lutArray = [pdata.calibrationVector.(lut)];
assert( length( lutArray ) / NrgCal == Ncal, ['The number of calibration ' ...
    'values for the desired LUT is not identical for all calibration ' ...
    'vectors. This function cannot proceed.'] );
lutArray = reshape( lutArray, [NrgCal Ncal] );

% construct the image LUT matrix, if necesary
if ~strcmp( imgLut, 'original' )
    imgLutArray = [pdata.calibrationVector.(imgLut)];
    assert( length( imgLutArray ) / NrgCal == Ncal, ['The number of ' ...
        'calibration values for the image LUT is not identical for all ' ...
        'calibration vectors. This function cannot proceed.'] );
    imgLutArray = reshape( imgLutArray, [NrgCal Ncal] );
end

% calculate relative azimuth time for calibration vectors
calTime = reshape( [pdata.calibrationVector.azimuthTime], [6, Ncal] ).';
calTimeRel = etime( calTime, T0 );
clear pdata calTime T0;


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% interpolate LUT data and apply it to input image

% create a mesh grid of calibration range pixels and azimuth times
[calAzGrid, calRgGrid] = meshgrid( calTimeRel, calRg );
clear calTimeRel calRg;

% create a mesh grid of interpolation range pixels and azimuth times
[imgAz, imgRg] = meshgrid( reshape( trel.', [1, Naz * Nb] ), (0:Nrg-1) );
clear trel;

% perform a bilinear interpolation of the desired LUT matrix
lutArrayInt = interp2( calAzGrid, calRgGrid, lutArray, imgAz, imgRg, ...
    'linear', NaN );
lutArrayInt = reshape( lutArrayInt, [Nrg, Naz, Nb] );
clear lutArray;

% perform a bilinear interpolation of the image LUT matrix, if necessary
if ~strcmp( imgLut, 'original' )
    imgLutArrayInt = interp2( calAzGrid, calRgGrid, imgLutArray, imgAz, ...
        imgRg, 'linear', NaN );
    imgLutArrayInt = reshape( imgLutArrayInt, [Nrg, Naz, Nb] );
    clear imgLutArray;
end
clear calAzGrid calRgGrid imgAz imgRg;

% display progress information, if requested
if verbose
    fprintf( 'Done!\n' );
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

% display progress information, if requested
if verbose
    fprintf( 'Done!\n' );
end

