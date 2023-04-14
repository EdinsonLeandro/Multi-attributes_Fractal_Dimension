function [dt,nsamples,ntraces] = checkSegyData(filenames)
% Check if the following data are the same for all input seismic volumes:
% Sampling interval.
% Number of samples.
% Number of traces in the seismic volume.
% Coordinates of the first and last trace.
% 
% Otherwise, the execution must be interrupted.
% "filenames" Matrix is a cell array that contains the names of the input files.

% Description of variables:
% segyHeader:       Segy Header of seismic volumes.
% segyTraceHeader:  Segy Trace Header of seismic traces.
% params:           Matrix containing the parameters to compare.
% checkSampleTrace: Check if the sampling interval, number of samples
%                   and number of traces are the same.
% checkCoordinate:  Check if the coordinates are the same.

% Link of interest:
% https://stackoverflow.com/questions/27535339/how-to-check-if-all-rows-of-a-matrix-are-equal

nFiles=size(filenames,2);

%% Check sampling interval, number of samples and number of traces

% Get parameters
params = zeros(nFiles,3);
for iFile = 1:nFiles
    [segyHeader] = ReadSegyHeader(filenames{iFile});
    params(iFile,:) = [segyHeader.dt/1000 ; segyHeader.ns ; segyHeader.ntraces];
end

% nnz() returns the number of nonzero matrix elements.
if nnz( diff(params, 1) ) ~= 0
    % There is a difference not equal to zero.
    error('Error in SegyHeader of seismic volumes');
end

%% Check coordinates

% Coordinates of the first and last trace.
cdpXY = zeros(nFiles,2,2);
for iFile = 1:nFiles
    % 'traces': Read specific time trace number
    % 'SkipData': Read only the header values (Data will return empty)
    [~ , segyTraceHeader] = ReadSegy(filenames{iFile}, 'traces', [1, params(1,3)], 'SkipData', 1);
    cdpXY(iFile,:,1) = [segyTraceHeader(1).cdpX , segyTraceHeader(1).cdpY];
    cdpXY(iFile,:,2) = [segyTraceHeader(2).cdpX , segyTraceHeader(2).cdpY];
end

% Check if the coordinates of the first and the last trace are equal
checkCoordinate = zeros(1, 2);
for h=1:2
    % Difference between consecutive elements across columns. "diff(.,1)"
    % Detect if all such differences are zeros. "nnz(..)==0"
    checkCoordinate(1,h) = nnz( diff(cdpXY(:,:,h), 1) );
end

% any() True if any element of a vector is a nonzero number.
if any(checkCoordinate(:))
    disp('Error in seismic volume coordinates');
end

% Return variables
dt = params(1,1);
nsamples = params(1,2);
ntraces = params(1,3);