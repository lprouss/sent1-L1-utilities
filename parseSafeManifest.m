function maniHead = parseSafeManifest( dsetPath )
% Extract information from the SAFE manifest file of a Sentinel-1 dataset.
%
% Inputs:
%   - dsetPath: path to the Sentinel-1 dataset directory
%
% Outputs:
%   - maniHead: structure containing the extracted information
%
% Required functions (not part of MATLAB):
%   - getXMLnode (see https://github.com/lprouss/xmlExtract)
%
% Additional information:
%   Although the SAFE format may be used for many ESA platforms, this function
%   only supports Level 1 Sentinel-1 products (SLC and GRD).
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: May 2014
% Updated: October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% try to read the input manifest file (XML tree)
maniFile = fullfile( dsetPath, 'manifest.safe' ); % manifest file
try
    % success
    xdoc = xmlread( maniFile );
catch err
    % fail, return exception
    rethrow( err );
end

%% check if the input file is indeed a manifest file
xroot = xdoc.getDocumentElement; % root node of the tree
if ~strcmp( xroot.getNodeName, 'xfdu:XFDU' )
    % if not, return error
    error( 'The input file is not a SAFE manifest file.' );
end


%% read metadataObject section
meta = xroot.getElementsByTagName('metadataSection').item(0);

%% get satellite platform name
platform = lower( getXMLnode( 'safe:familyName', meta, 'str' ) );
% capitalize first letter
platform = regexprep( platform, '(\<\w)', '${upper($1)}' );
% make sure the platform is indeed Sentinel-1
assert( strcmpi( platform, 'sentinel-1' ), ['The platform in the manifest ' ...
    'file is not supported: %s'], platform );
% add platform number to the name
platform = strcat( platform, getXMLnode( 'safe:number', meta, 'str' ) );

%% get orbit pass
pass = lower( getXMLnode( 's1:pass', meta, 'str' ) );

%% get product type
prodType = getXMLnode( 's1sarl1:productType', meta, 'str' );
% make sure the product is SLC or GRD
prodTest = any( strcmpi( prodType, {'SLC', 'GRD'} ) ); % test for product type
assert( prodTest, ['The product type in the manifest file is not ' ...
    'supported: %s'], prodType );
clear prodTest;

%% get acquisition mode
acqMode = getXMLnode( 's1sarl1:mode', meta, 'str' );

%% get number of subswaths and IDs
% node for swath(s)
node = meta.getElementsByTagName('s1sarl1:swath');
numSwath = node.getLength; % number of subswaths
swathID = cell( 1, numSwath ); % allocate memory for subswaths IDs
% get IDs for all subswaths
for cnt = 1:numSwath
    swathID{cnt} = char( node.item(cnt-1).getTextContent );
end

%% get polarizations
% node for polarization(s)
node = meta.getElementsByTagName('s1sarl1:transmitterReceiverPolarisation');
numPol = node.getLength; % number of polarizations
pol = cell( 1, numPol ); % allocate memory for polarizations
% get polarization name(s)
for cnt = 1:numPol
    pol{cnt} = char( node.item(cnt-1).getTextContent );
end

%% number of data channels
numChan = length( pol ) * length( swathID );

%% get acquisition start and stop times
dateStrFmt = 'yyyy-mm-ddTHH:MM:SS.FFFFFF';
acqT0 = getXMLnode( 'safe:startTime', meta, 'dateStr', dateStrFmt );
acqTE = getXMLnode( 'safe:stopTime', meta, 'dateStr', dateStrFmt );

%% read dataObject section
data = xroot.getElementsByTagName('dataObjectSection').item(0);

%% allocate memory for data and annotation files and initialize counters
dataFile = cell( 1, numChan ); cntD = 1;
prodAnn = cell( 1, numChan ); cntP = 1;
calAnn = cell( 1, numChan ); cntC = 1;
noiseAnn = cell( 1, numChan ); cntN = 1;

%% get data and annotations files names (attributes)
for cnt = 1:data.getLength
    % only process current item if it is not empty
    if data.item(cnt-1).hasAttributes
        % repID attribute name
        repID = char( data.item(cnt-1).getAttribute('repID') );
        % get fileLocation node
        loc = data.item(cnt-1).getElementsByTagName('fileLocation').item(0);
        % determine the type of file
        if strcmp( repID, 's1Level1MeasurementSchema' )
            % measurement dataset
            str = char( loc.getAttribute('href') );
            dataFile{cntD} = str(3:end);
            cntD = cntD + 1;
        elseif strcmp( repID, 's1Level1ProductSchema' )
            % product annotation
            str = char( loc.getAttribute('href') );
            prodAnn{cntP} = str(3:end);
            cntP = cntP + 1;
        elseif strcmp( repID, 's1Level1CalibrationSchema' )
            % calibration annotation
            str = char( loc.getAttribute('href') );
            calAnn{cntC} = str(3:end);
            cntC = cntC + 1;
        elseif strcmp( repID, 's1Level1NoiseSchema' )
            % noise annotation
            str = char( loc.getAttribute('href') );
            noiseAnn{cntN} = str(3:end);
            cntN = cntN + 1;
        end
    end
end


%% construct the manifest header structure
maniHead.rootDir = fileparts( maniFile ); % root directory
maniHead.platform = platform; % satellite platform
maniHead.pass = pass; % orbit pass
maniHead.prodType = prodType; % product type
maniHead.acqMode = acqMode; % acquisition mode
maniHead.swathID = swathID; % subswaths IDs
maniHead.polarization = pol; % polarizations
maniHead.acqStartTime = acqT0;
maniHead.acqStopTime = acqTE;
maniHead.dataFile = dataFile; % data files names
maniHead.prodAnn = prodAnn; % product annotation files names
maniHead.calAnn = calAnn; % calibration annotation files names
maniHead.noiseAnn = noiseAnn; % noise annotation files names

