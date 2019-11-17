function data = readSent1Data( dataFile, varargin )
% Read data in a Sentinel-1 TIFF (measurement) file.
%
% Inputs:
%   - dataFile: full path and name of the Sentinel-1 image file to read.
%   - roi (optional): structure defining the region of interest (ROI) to read
%       in the data file, see Additional information.
%   - verbose (optional): logical flag to display progress information, false
%       by default.
%
% Outputs:
%   - data: data array read in the input Sentinel-1 data file, with dimensions
%       [Nrg, Naz] where 'Nrg' and 'Naz' are the number of range samples and
%       azimuth lines read in the file, respectively.
%
% Required functions (toolboxes and/or user-defined):
%   - getTIFFinfo
%
% Additional information:
%   The ROI to read in data is defined using the following parameters in the
%   input 'roi' structure:
%   - 'ridx': indexes of range samples to read.
%   - 'aidx': indexes of azimuth lines to read.
%   Indexes should start at 1. Default value is 0 for all parameters, i.e. all
%   data is read in that dimension. When 'roi' is omitted, the entire data file
%   is read.
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: May 2014
% Updated: November 2017, November 2019 (simplified ROI structure, other minor
%   improvements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% initialization

% validate the number of inputs
narginchk( 1, 3 );

% open the input TIFF data file for reading and get information to read it
[fileInfo, fid] = getTIFFinfo( dataFile );

% assign variable to the provided optional input(s)
roi = struct(); % empty ROI structure
verbose = false; % no verbose
for ni = 1:nargin-1
    % assign current optional input to correct variable using its type
    if isstruct( varargin{ni} )
        % structure: reading options
        roi = varargin{ni};
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
%%% verify provided ROI structure or create a default one

% number of range samples and azimuth lines in the data file
Nrg = fileInfo.numPixels;
Naz = fileInfo.numLines;

% indexes of range samples and azimuth lines to read for ROI
if isfield( roi, 'ridx' )
    ridx = unique( roi.ridx(:) );
    validateattributes( ridx, {'numeric'}, {'integer', 'nonnegative', ...
        '<=', Nrg}, '', 'roi.ridx' );
    if any( ridx ) == 0
        ridx = ( 1:Nrg );
    end
else
    % read all range samples in each swath
    ridx = ( 1:Nrg );
end
if isfield( roi, 'aidx' )
    aidx = unique( roi.aidx(:) ).';
    validateattributes( aidx, {'numeric'}, {'integer', 'nonnegative', ...
        '<=', Naz}, '', 'roi.aidx' );
    if any( aidx ) == 0
        aidx = ( 1:Naz );
    end
else
    % read all azimuth lines in each swath
    aidx = ( 1:Naz );
end
NrgRoi = numel( ridx );
NazRoi = numel( aidx );


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% read the data file

% display progress information, if requested
if verbose
    fprintf( 'Reading data from the file' );
end

% calculate bytes offsets for the ROI
roiBytesOffset = fileInfo.linesBytesOffset(aidx) + ( ridx(1) - 1 ) * ...
    fileInfo.bytesPerSample * ( 1 + fileInfo.complexFlag );

% dimensions of the data array to read
dataDim = [NrgRoi * ( 1 + fileInfo.complexFlag ), NazRoi];

% calculate number of bytes to skip after reading each line
NbytesLine = NrgRoi * ( 1 + fileInfo.complexFlag ) * ...
    fileInfo.bytesPerSample; % number of bytes per line
NbytesSkip = roiBytesOffset(2:end) - ( roiBytesOffset(1:end-1) + NbytesLine );
clear NbytesLine;

% data type (save as single)
dtype = fileInfo.dataType;
if ~strcmpi( fileInfo.dataType, 'single' )
    dtype = [dtype, '=>single'];
end

% read data in file
if numel( unique( NbytesSkip ) ) == 1
    % number of bytes to skip after each line is constant
    NbytesSkip = unique( NbytesSkip );

    % display progress information, if requested
    if verbose
        fprintf( ' (all data read at once)... ' );
    end

    % set a pointer at the beginning of the ROI in the file
    try
        % try to set pointer position
        fseek( fid, roiBytesOffset(1), -1 );
    catch
        % failed setting pointer position, return error
        error( ['Error setting the pointer at the right position in data ' ...
            'file: %s.'], dataFile );
    end

    % add the number of samples to read to the data type string if necessary
    if NbytesSkip > 0
        % number of bytes to skip is not zero, modify string
        dtype = sprintf( '%d*%s=>single', NrgRoi * ( 1 + ...
            fileInfo.complexFlag ), fileInfo.dataType );
    end

    % read all samples in the ROI simultaneously and convert to single data type
    data = fread( fid, dataDim, dtype, NbytesSkip );
    clear dtype;
else
    % number of bytes to skip after each line is not constant
    if verbose
        fprintf( ':\n' );
    end

    % initialize an array for the ROI data
    data = zeros( dataDim, 'single' );

    % read data for the ROI line-by-line
    for nl = 1:NazRoi
        % display progress information, if requested
        if verbose
            perc = round( ( nl - 1 ) / NazRoi * 100 );
            fprintf( '  %d%% complete\r', perc );
            clear perc;
        end

        % set a pointer at the beginning of the current line in the file
        try
            % try to set pointer position
            fseek( fid, roiBytesOffset(nl), -1 );
        catch
            % failed setting pointer position, return error
            error( ['Error setting the pointer at the right position in data ' ...
                'file: %s.'], dataFile );
        end

        % read current line in the file
        data(:,nl) = fread( fid, NrgRoi * ( 1 + fileInfo.complexFlag ), ...
            dtype );
    end
    if verbose
        fprintf( '\n' );
    end
end
clear NbytesSkip roiBytesOffset;

% close the open data file
fclose( fid );

% if necessary, convert samples to complex data
if fileInfo.complexFlag
    %data = data(1:2:end,:) + 1j * data(2:2:end,:);
    data = complex( data(1:2:end,:), data(2:2:end,:) );
end
clear fileInfo;

% display progress information, if requested
if verbose
    fprintf( 'Done!\n' );
end

