function [img, varargout] = mergeSent1TopsBursts( data, tvecB );
% Stitch (merge) burst images for in a Sentinel-1 TOPS subswath.
%
% Inputs:
%   - data: data array containing the burst images for the TOPS subswath. It
%       should have the following dimensions: [Nrg Nlb Nb], where 'Nrg' is the
%       number of range pixels, 'Nlb' is the number of azimuth pixels in each
%       burst image and 'Nb' is the number of bursts.
%   - tvecB: matrix of azimuth time vectors for the bursts, with
%       dimensions [Nb Nlb].
%
% Outputs:
%   - img: stitched TOPS image, with dimensions [Nrg Naz], where 'Naz' is the
%       number of azimuth pixels after bursts stitching.
%   - imgPars (optional): information structure for bursts merging.
%
% Required functions (not part of MATLAB): none
%
% Author: Louis-Philippe Rousseau (Université Laval)
% Created: January 2018
% Updated: May 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% validate dimensions of the input data array
assert( ndims( data ) == 3, ['The input data array does not have the ' ...
    'expected number of dimensions. It should be a 3D array: ' ...
    '[Nrg x Nlb x Nb].'] );

% dimensions of data array
dataDim = size( data );
Nrg = dataDim(1);
Nlb = dataDim(2);
Nb = dataDim(3);
clear dataDim;

% total number of azimuth lines in the input data array
totLine = Nlb * Nb;

% validate dimensions of azimuth time vector array
assert( all( size( tvecB ) == [Nb, Nlb] ), ['Dimensions of the provided matrix ' ...
    'of azimuth time vectors for the bursts do not match the of the input data ' ...
    'array. It should be: [Nb x Nlb].'] );

% determine the azimuth line time interval in the burst images
lti = diff( tvecB, 1, 2 );
lti = unique( lti(:) ); % unique values
assert( any( lti - mean( lti ) ) > 1e-6, ['The azimuth line time interval ' ...
    'seems to vary in the input azimuth time matrix. This function cannot ' ...
    'proceed.'] );
if length( lti ) > 1
    % if more than one valid LTI value, use average
    lti = mean( lti );
end

% bursts duration [s]
%Tburst = lti * Nlb;
%Tburst = lti * ( Nlb - 1 );

% number of azimuth lines in the merged image
azIdxImg0 = round( tvecB(1) / lti );
azIdxImgE = round( tvecB(end) / lti );
Naz = azIdxImgE - azIdxImg0 + 1;

% azimuth time vector for the merged image
tvecImg = ( azIdxImg0 + ( 0:Naz-1 ) ) * lti;

% make sure the azimuth time origin is zero for next calculations
tvecB = tvecB - min( tvecB(:) );

% relative start and stop times of bursts
%tb0Rel = tvecB(:,1).';
%tbEndRel = tvecB(:,end).';

% overlap between successive bursts
%Tover = tbEndRel(1:end-1) - tb0Rel(2:end); % temporal overlap
Tover = tvecB(1:end-1,end) - tvecB(2:end,1);
Tover = [0; Tover; 0]; % add a '0' overlap for first and last bursts
Nover = Tover / lti; % number of overlap lines
%Nover = round( Tover / lti ); % number of overlap lines
%NoverHalf = Nover / 2; % half the number of overlap lines (real value)
%clear tb0Rel tbEndRel Tover;

% start azimuth indexes for subset of burst images
azIdxB0 = round( Nover(1:end-1) / 2 );

% start azimuth indexes for burst images in merged array
%azIdxImg0 = round( tvecB(:,1) / lti ) + floor( Nover(1:end-1) / 2 );
azIdxImg0 = round( tvecB(:,1) / lti ) + azIdxB0;
%azIdxImg0 = round( ( tvecB(:,1) + Tover(1:end-1) / 2 ) / lti );

% number of azimuth samples in subsets of burst images
NazSubB = [ diff( azIdxImg0 ); Nlb - azIdxB0(end)];

% number of azimuth lines in stitched image
%Naz = totLine - sum( Nover );

% start and stop azimuth indexes to subset data in each burst
%b0Idx = floor( NoverHalf(1:end-1) ) + 1;
%bEndIdx = Nlb - ceil( NoverHalf(2:end) );
%clear NoverHalf Nover;
%NazSubB = Nlb - azIdxB0 - ceil( Nover(2:end) / 2 );
%azIdxB = azIdxB0 + ( 1:NazSubB );

% start and stop azimuth indexes for each burst in the stitched image
%imgEndIdx = cumsum( bEndIdx - b0Idx + 1 );
%img0Idx = [0 imgEndIdx(1:end-1)] + 1;
%azIdxImg = azIdxImg0 + ( 1:NazSubB );

% initialize array for the stitched image
img = zeros( Nrg, Naz, class( data ) );

% initialize information structure for bursts merging
imgPars.tazImg = zeros( 1, Naz );
imgPars.burstIdx = cell( 1, Nb );
imgPars.imgIdx = cell( 1, Nb );

% merge the bursts
for nb = 1:Nb
    % azimuth indexes for the original data array
    %dataIdx = ( b0Idx(nb):bEndIdx(nb) );
    dataIdx = azIdxB0(nb,:) + ( 1:NazSubB(nb) );
    imgPars.burstIdx{nb} = dataIdx;

    % azimuth indexes for the image array
    %imgIdx = ( img0Idx(nb):imgEndIdx(nb) );
    imgIdx = azIdxImg0(nb,:) + ( 1:NazSubB(nb) );
    imgPars.imgIdx{nb} = imgIdx;

    % copy burst data from the original array into the image array
    img(:,imgIdx) = data(:,dataIdx,nb);
    imgPars.tazImg(imgIdx) = tvecB(nb,dataIdx);
    clear imgIdx dataIdx;
end
%clear data b0Idx bEndIdx img0Idx imgEndIdx;

% validate azimuth time vector for merged image
assert( all( abs( imgPars.tazImg - tvecImg ) < 1e-6 ), ['Problem with bursts ' ...
    'merging: mismatch between expected and real azimuth time vectors after ' ...
    'merging.'] );

% return parameters structure for processed image, if requested
if nargout == 2
    varargout{1} = imgPars;
end

