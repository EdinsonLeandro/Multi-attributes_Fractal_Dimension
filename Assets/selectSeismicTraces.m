function seismicTraces = selectSeismicTraces(N, surface, xyNumTraces)
% Select X-Y coordinates and the number of seismic traces where surface
% data was interpolated.
% This will be the traces to calculate Fractal Dimension.
% Input:
% N:           Number of seismic traces in segy data.
% surface:     Interpolated surface data. X-Y-Time data.
% xyNumTraces: Matrix with X-Y coordinates of seismic traces of all
%              seismic volume and number of traces.

% Description of the variables:
% seismicTraces. Matrix with X-Y coordinates and the number of sekected seismic
%                traces. The size will be lower than or equal to the number
%                of traces in seismic volume.
% Assumption: The input surface must be interpolated in all area of seismic data.

% Initizalice "seismicTraces". Empty cells will be removed.
seismicTraces = zeros(N, 3);
kTrace = 1;

% Instead of looking for coordinate by coordinate in surface data,
% we use a more efficient method. For each unique X-coordinate,
% find Y-coordinates and save data.

% Unique values for X coordinates in surface data.
xCoordSurface = unique(surface(:, 1), 'stable');

% Iterate over each X-coordinate in surface data
for iXCoord = 1 : size(xCoordSurface, 1)

    % Find indexes in "xyNumTraces" where the X-coordinate of seismic trace
    % are equal to X-coordinate in surface data.
    indexesXcoord = find(xyNumTraces(:, 1) == xCoordSurface(iXCoord));
    
    % If "indexesXcoord" is empty, go for the next X-coordinate in surface data.
    if ~isempty(indexesXcoord)
        
        % Both X-coordinates match. Select all Y-coordinates of seismic traces.
        yCoordTrace = xyNumTraces(indexesXcoord, 2);
        
        % Select Y-coordinates of surface data
        yCoordSurface = surface(surface(:, 1) == xCoordSurface(iXCoord), 2);
        
        % Iterate over each Y-coordinate in surface data only in X-coordinate
        % previously selected
        for jYCoord = 1:size(yCoordSurface)
            % First index where Y-coordinates of seismic traces is equal to
            % Y-coordinate in surface data selected
            firstIndex = find(yCoordTrace == yCoordSurface(jYCoord), 1);
            
            % Save X-coordinate, Y-coordinate and number of trace
            seismicTraces(kTrace, :) = xyNumTraces(indexesXcoord(1) + firstIndex - 1, :);
            
            kTrace = kTrace + 1;
            
        end % End for

    end % End if
    
end % End for

% Remove zero rows
seismicTraces( ~any(seismicTraces, 2), : ) = [];