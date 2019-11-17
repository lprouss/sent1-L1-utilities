function [tiffInfo, varargout] = getTIFFinfo( inpFile )
% Extract parameters from the specified TIFF file, perform basic checks and
% return extracted information and a pointer to the open TIFF file, if
% requested.
%
% Inputs:
%   - inpFile: TIFF file name.
%
% Outputs:
%   - tiffInfo: structure containing extracted image parameters.
%   - fid (optional): pointer to the open TIFF file.
%
% Required functions (not part of MATLAB):
%   - getTIFFfield (internal)
%
% Additional information:
%   This function reads information from the TIFF tags in the file and uses them
%   to extract parameters. It is assumed that the TIFF file follows the TIFFv6.0
%   specification (see
%   http://partners.adobe.com/public/developer/tiff/index.html for details).
%
%   This function supports RADARSAT-2 images in the BigTIFF format, Sentinel-1
%   imagery, and RADARSAT-2's MODEX-1/2 (experimental mode) data processed
%   using the DRDC/MDA GMTI processor (phase 3) at DRDC-Ottawa.
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: May 2014
% Updated: June 2014, October 2017, November 2017, January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TODO: add support for GeoTIFF tags?

% open input file for reading and check for validity
fid = fopen( inpFile, 'r', 'ieee-le' ); % first open as little-endian
assert( fid >= 0, 'Failed to open the input file for reading: %s', inpFile );

% read the first two bytes (0-1) in the file: byte order
order = fread( fid, 2, 'uint8=>char' ); % read and transform to char
if strcmp( order', 'MM' )
    % the two characters are 'MM' (big-endian), re-open file as such
    fclose( fid ); % first close the file
    fid = fopen( inpFile, 'r', 'ieee-be' ); % open as big-endian
    byteOrder = 'be'; % set file byte order as big-endian
elseif strcmp( order', 'II' )
    % the two characters are 'II' (little-endian), nothing to do
    byteOrder = 'le'; % set file byte order as little-endian
else
    % otherwise, close the file and return error
    fclose( fid );
    error( 'Invalid byte ordering in the input file.' );
end

% get the byte order of the system running the function
[~, ~, endian] = computer;

% set a flag to change the data byte ordering to that of the system
boFlag = ( (endian == 'L' && order(1) == 'M') || ...
    (endian == 'B' && order(1) == 'I') );
clear order endian;

% read the next two bytes (2-3) in the file: file type
fileType = fread( fid, 1, 'int16' ); % read as int16

% read the image file directory (IFD), according to the file type
switch fileType
    case 42
        % the file is in the standard TIFF format

        % set data type for numbers and offsets
        offsetType = 'uint32';

        % get offset from header (bytes 4-7)
        ifdOffset = fread( fid, 1, 'uint32' ); % read as uint32
    case 43
        % the file is in the BigTIFF format

        % get the bytes size of offsets (bytes 4-5, should be 8)
        offsetBytes = fread( fid, 1, 'int16' ); % read as int16
        assert( offsetBytes == 8, ['The bytes size of offsets is %d, but 8 ' ...
            'is expected for BigTIFF files.'], offsetBytes );
        clear offsetBytes;

        % make sure the next field is 0 (bytes 6-7)
        val = fread( fid, 1, 'int16' ); % read as int16
        assert( val == 0, ['The field after the bytes size of offsets should ' ...
            'be 0 in BigTIFF files, but is %d instead.'], val );
        clear val;

        % set data type for numbers and offsets
        offsetType = 'uint64';

        % get offset from header (bytes 8-15)
        ifdOffset = fread( fid, 1, 'uint64' ); % read as uint64
    otherwise
        % the file is not in TIFF format, close it and return error
        fclose( fid );
        error( 'The input file is not in TIFF format.' );
end

% move to the position of the IFD in the file
val = fseek( fid, ifdOffset, -1 );
assert( val >= 0, ['Error positioning the pointer at the image file ' ...
    'directory (IFD) in the TIFF file.'] );
clear val;

% read the number of fields in the IFD
if fileType == 42
    % standard TIFF file, read the first two bytes
    numFields = fread( fid, 1, 'uint16' ); % read as uint16
elseif fileType == 43
    % BigTIFF file, read the first eight bytes
    numFields = fread( fid, 1, 'uint64' ); % read as uint64
end

% allocate the memory for IFD's fields data
% four entries per field: tag, value(s) type, value(s) count, value(s) offset
ifdfield(numFields).tag = [];
ifdfield(numFields).type = [];
ifdfield(numFields).count = [];
ifdfield(numFields).offset = [];

% read the IFD's fields
for cnt = 1:numFields
    % the first two fields are read as uint16
    ifdField(cnt).tag = fread( fid, 1, 'uint16' );
    ifdField(cnt).type = fread( fid, 1, 'uint16' );

    % the last two fields are read as uint32 (TIFF) or uint64 (BigTIFF)
    ifdField(cnt).count = fread( fid, 1, offsetType );
    ifdField(cnt).offset = fread( fid, 1, offsetType );
end
clear numFields;

% get data for additional IFD fields
% number of pixels per line (tag 256)
numPixels = getTIFFfield( fid, ifdField, 256, offsetType, boFlag );
% number of lines (tag 257)
numLines = getTIFFfield( fid, ifdField, 257, offsetType, boFlag );
% number of bits per sample (tag 258)
numBitsSamp = getTIFFfield( fid, ifdField, 258, offsetType, boFlag );
% bytes offset for each data strip (tag 273)
offsetStrips = getTIFFfield( fid, ifdField, 273, offsetType, boFlag );
% number of samples per pixel (tag 277)
numSampPixel = getTIFFfield( fid, ifdField, 277, offsetType, boFlag );
% number of rows per strip (tag 278)
numRowsStrips = getTIFFfield( fid, ifdField, 278, offsetType, boFlag );
% number of bytes in each data strip (tag 279)
%numBytesStrips = getTIFFfield( fid, ifdField, 279, offsetType, boFlag );
% sample format (tag 339)
%SampFormat = getTIFFfield( fid, ifdField, 339, offsetType, boFlag );

% set data type
numBytesSamp = numBitsSamp / 8; % number of bytes per sample
switch numBytesSamp
    case 1
        % 1 byte per sample
        dataType = 'uint8';
        cplxFlag = 0;
    case 2
        % 2 bytes per sample
        if numSampPixel == 1
            % 1 sample per pixel (magnitude or power data)
            dataType = 'uint16';
            cplxFlag = 0;
        elseif numSampPixel == 2
            % 2 samples per pixel: RADARSAT-2's SLC/raw (complex) data
            dataType = 'int16';
            cplxFlag = 1;
        else
            % unsupported number of samples per pixel, return error
            error( 'Unsupported number of samples per pixel, %d.', ...
                numSampPixel );
        end
    case 4
        % 4 bytes per sample: Sentinel-1 SLC (complex) data
        numBytesSamp = 2;
        dataType = 'int16';
        cplxFlag = 1;
    case 8
        % 8 bytes per sample: RADARSAT-2 MODEX (complex) data
        numBytesSamp = 4;
        %dataType = 'int32';
        dataType = 'single';
        cplxFlag = 1;
    otherwise
        % unsupported number of bytes per sample, return error
        error( 'Unsupported number of bytes per sample, %d.', numBytesSamp );
end
clear numSampPixel numBitsSamp;

% calculate the bytes offset for each row in a data strip
offsetRows = (0:numRowsStrips-1)' * numBytesSamp * (1 + cplxFlag) * numPixels;
%clear numRowsStrips;

% calculate the bytes offset for each data line
linesBytesOffset = bsxfun( @plus, offsetRows, offsetStrips(:)' );
linesBytesOffset = linesBytesOffset(:)'; % reshape as a line
clear offsetRows offsetStrips;


% create the info structure
tiffInfo.endian = byteOrder;
%tiffInfo.numPixels = single( numPixels );
%tiffInfo.numLines = single( numLines );
tiffInfo.numPixels = numPixels;
tiffInfo.numLines = numLines;
tiffInfo.dataType = dataType;
tiffInfo.complexFlag = cplxFlag;
tiffInfo.bytesPerSample = numBytesSamp;
tiffInfo.linesBytesOffset = linesBytesOffset;
clear byteOrder numPixels numLines cplxFlag dataType numButesSamp ...
    linesBytesOffset;


% return pointer to open TIFF file or close it
if nargout == 2
    % pointer requested, return it to caller
    varargout{1} = fid;
else
    % pointer not requested, close the file
    fclose( fid );
end

%%% end of 'getTIFFinfo'
end


%%% beginning of internal function 'getTIFFfield'
function value = getTIFFfield( fid, ifdData, tag, offsetType, boFlag )
% Get the value or array of values for the TIFF field with the specified
% tag, using the information in the image file directory (IFD).
%
% Inputs:
%   - fid: pointer to the open TIFF file.
%   - ifdData: structure containing information about the IFD.
%   - tag: desired field tag (decimal).
%   - offsetType: data type for offsets.
%   - boFlag: flag to change the data byte ordering.
%
% Outputs:
%   - value: extracted value(s) for the TIFF field.
%
% Required functions (not part of MATLAB): none
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: May 2014
% Updated: June 2014, October 2017, Novermber 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the requested tag in the IFD data structure
idx = find( [ifdData.tag] == tag ); % index of the tag
assert( ~isempty( idx ), 'Tag %s is not present in the IFD.', num2str( tag ) );

% check the number of values (count) for the field
if ifdData(idx).count == 1
    % only one value for the field, the offset is the value
    value = ifdData(idx).offset;

    % set the byte size of value (4 or 8)
    if strcmp( offsetType, 'uint32' )
        % 4 bytes offsets (TIFF)
        bsize = 4;
    elseif strcmp( offsetType, 'uint64' )
        % 8 bytes offsets (BigTIFF)
        bsize = 8;
    else
        % invalid offset type, return error
        error( 'Invalid bytes size for offsets, %s.', offsetType );
    end

    % shift bits, if the flag for byte order change is set
    if boFlag
        switch ifdData(idx).type
            case {1, 2}
                % data type is uint8 or ascii, shift 3 or 7 bytes to the left
                value = bitshift( value, -(bsize-1)*8 );
            case 3
                % data type is uint16, shift 2 or 6 bytes to the left
                value = bitshift( value, -(bsize-2)*8 );
            case 4
                % data type is uint32, shift 0 or 4 bytes to the left
                value = bitshift( value, -(bsize-4)*8 );
            otherwise
                % data type is uint64 (BigTIFF only), do not shift
        end
    end
    clear bsize;
else
    % array of values for the field

    % use offset to move to the right position in the file
    val = fseek( fid, ifdData(idx).offset, -1 );
    assert( val >= 0, ['Error setting the pointer at the right position in ' ...
        'the TIFF file.'] );
    clear val;

    % read the values according to count and data type
    switch ifdData(idx).type
        case 1
            % data type is uint8
            value = fread( fid, ifdData(idx).count, 'uint8' );
        case 2
            % data type is ascii
            value = fread( fid, ifdData(idx).count, 'int8' );
        case 3
            % data type is uint16
            value = fread( fid, ifdData(idx).count, 'uint16' );
        case 4
            % data type is uint32
            value = fread( fid, ifdData(idx).count, 'uint32' );
        case 16
            % data type is uint64 (BigTIFF only)
            value = fread( fid, ifdData(idx).count, 'uint64' );
        otherwise
            % data type is not supported, return error
            error( 'Unsupported data type %d', ifdData(idx).type );
    end
end
clear idx;

%%% end of internal function 'getTIFFfield'
end

