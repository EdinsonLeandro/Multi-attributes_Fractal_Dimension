function coordinates = traceNumbers(file, n)
% Linear sequence of seismic traces.
% Input:
% file: Segy file.
% n: Number of seismic traces in segy file.

% Open and read all trace headers from segy data. In this line is where
% the function takes de longest time.
[~ , allTraceHeaders] = ReadSegy(file, 'SkipData',1);

% Coordinates of seismic traces and number of traces
coordinates = zeros(n, 3);

for iTrace = 1:n
    % Divide each coordinate by 100 to bring them to metric scale
    coordinates(iTrace, :) = [round(allTraceHeaders(iTrace).SourceX / 100), ...
                              round(allTraceHeaders(iTrace).SourceY / 100), ...
                              allTraceHeaders(iTrace).TraceNumber];
end

% Get unique X coordinates
uniqueCoord = unique(coordinates(:, 1), 'stable');

% Number of traces by seismic line
totalTraces = zeros(size(uniqueCoord));
for index = 1:size(uniqueCoord,1)
    % Find indexes where X coordinate is equal to each value
    condition = (coordinates(:, 1) == uniqueCoord(index, 1));
    
    % Save the total number of traces by Inline
    totalTraces(index) = size(coordinates(condition), 1);

    % Update trace number with linear sequence
    coordinates(condition, 3) = coordinates(condition, 3) + (totalTraces(1) * (index-1));
end

if ~all(totalTraces == totalTraces(1))
    error('The number of traces by seismic Inline is not the same')   
end
