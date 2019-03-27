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
% Required functions (not part of MATLAB):
%   - getTIFFinfo
%
% Additional information:
%   The ROI to read in data is defined using the following parameters in the
%   input 'roi' structure:
%   - roi.firstSample: index of the first range sample in the desired ROI
%       (starts at 1);
%   - roi.numSample: number of range samples in the desired ROI.
%   - roi.firstLine: index of the first azimuth line in the desired ROI
%       (starts at 1);
%   - roi.numLine: number of azimuth lines in the desired ROI;
%   When 'roi' is omitted, the entire image is read from the data file.
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: May 2014
% Updated: November 2017
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

% field names and defaults for the read options structure
roiFields = {'firstLine', 1;
             'numLine', Naz;
             'firstSample', 1;
             'numSample', Nrg};

% make sure all fields are present in the options structure
fieldFlag = isfield( roi, roiFields(:,1) ); % flag provided fields
missingFields = roiFields(~fieldFlag,:); % list of missing fields
for cnt = 1:size( missingFields, 1 )
    % assing default values to missing fields
    roi.(missingFields{cnt,1}) = missingFields{cnt,2};
end
clear fieldFlag missingFields roiFields;

% validate or set azimuth lines window
if roi.firstLine <= 0 || roi.firstLine > Naz
    % invalid first line index, return warning and use default (1)
    warning( ['The provided index for the first azimuth line in the ROI is ' ...
        'invalid: %d. Using default (first sample).'], roi.firstLine );
    roi.firstLine = 1;
end
if roi.numLine < 1 || roi.firstLine - 1 + roi.numLine > Naz
    % invalid number of lines, return warning and use default
    NazRoi = Naz - roi.firstLine + 1; % number of lines to use
    warning( ['The provided number of azimuth lines in the ROI is invalid: ' ...
        '%d. Using %d lines: %d (total number of azimuth lines) - %d (first ' ...
        'line index) + 1.'], roi.numLine, NazRoi, Naz, roi.firstLine );
    roi.numLine = NazRoi; clear NazRoi;
end

% validate or set range samples window
if roi.firstSample <= 0 || roi.firstSample > Nrg
    % invalid first sample index, return warning and use default (1)
    warning( ['The provided index for the first range sample in the ROI is ' ...
        'invalid: %d. Using default (first sample).'], roi.firstSample );
    roi.firstSample = 1;
end
if roi.numSample < 1 || roi.firstSample - 1 + roi.numSample > Nrg
    % invalid number of samples, return warning and use default
    NrgRoi = Nrg - roi.firstSample + 1; % number of samples to use
    warning( ['The provided number of range samples in the ROI is invalid: ' ...
        '%d. Using %d samples: %d (total number of range samples) - %d ' ...
        '(first sample index) + 1.'], roi.numSample, NrgRoi, Nrg, ...
        roi.firstSample );
    roi.numSample = NrgRoi; clear NrgRoi;
end
clear Nrg Naz;


%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %%%
%%% read the data file

% display progress information, if requested
if verbose
    fprintf( 'Reading data from the file' );
end

% calculate bytes offsets for the ROI
roiLines = roi.firstLine + (0:roi.numLine-1); % list of lines in the ROI
roiBytesOffset = fileInfo.linesBytesOffset(roiLines) + ...
    (roi.firstSample - 1) * fileInfo.bytesPerSample * (1 + ...
    fileInfo.complexFlag);
clear roiLines;

% dimensions of the data array to read
dataDim = [roi.numSample * (1 + fileInfo.complexFlag), roi.numLine];

% calculate number of bytes to skip after reading each line
NbytesLine = roi.numSample * (1 + fileInfo.complexFlag) * ...
    fileInfo.bytesPerSample; % number of bytes per line
NbytesSkip = roiBytesOffset(2:end) - (roiBytesOffset(1:end-1) + NbytesLine);
clear NbytesLine;

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
        dtype = sprintf( '%d*%s', roi.numSample * (1 + ...
            fileInfo.complexFlag), fileInfo.dataType );
    else
        % number of bytes to skip is zeroo, do not change the string
        dtype = fileInfo.dataType;
    end

    % read all samples in the ROI simultaneously and convert to single data type
    data = single( fread( fid, dataDim, dtype, NbytesSkip ) );
    clear dtype;
else
    % number of bytes to skip after each line is not constant
    if verbose
        fprintf( ':\n' );
    end

    % initialize an array for the ROI data
    data = zeros( dataDim, fileInfo.dataType );

    % read data for the ROI line-by-line
    for nl = 1:roi.numLine
        % display progress information, if requested
        if verbose
            perc = round( (nl - 1) / roi.numLine * 100 );
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
        data(:,nl) = fread( fid, roi.numSample * (1 + fileInfo.complexFlag), ...
            fileInfo.dataType );
    end

    % convert data array to single data type
    data = single( data );
    if verbose
        fprintf( '\n' );
    end
end
clear NbytesSkip roiBytesOffset;

% close the open data file
fclose( fid );

% if necessary, convert samples to complex data
if fileInfo.complexFlag
    data = data(1:2:end,:) + 1j * data(2:2:end,:);
end
clear fileInfo;

% display progress information, if requested
if verbose
    fprintf( 'Done!\n' );
end

