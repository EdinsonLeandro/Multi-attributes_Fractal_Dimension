function filteredData = selectEqualCoord(N, reference, dataToFilter)
% Input:
% N:            Integer. Maximum number of samples.
% reference:    Matrix. Format X-Y-Value.
%               Use this data as reference to find the valid values.
% dataToFilter: Matrix. Format X-Y-Value.
%               Data to be filtered.
% Select X-Y coordinates and values where coordinates of "dataToFilter"
% are equal to the coordinates of the number of "reference".
%
% Ouput:
% filteredData. Matrix. Format X-Y-Value.
%               Data selected. The size will be lower than or equal to N.

% Initizalice "filteredData". Empty cells will be removed.
filteredData = zeros(N, 3);
kTrace = 1;

% Instead of looking for coordinate by coordinate in "reference" data,
% we use a more efficient method. For each unique X-coordinate,
% find Y-coordinates and save data.

% Unique values for X coordinates in "reference" data.
xCoordRef = unique(reference(:, 1), 'stable');

% Iterate over each X-coordinate in "reference" data.
for iXCoord = 1 : size(xCoordRef, 1)

    % Find indexes in "dataToFilter" where its X-coordinate are equal to
    % X-coordinate in "reference" data.
    indexesXcoord = find(dataToFilter(:, 1) == xCoordRef(iXCoord));
    
    % If "indexesXcoord" is empty, go for the next X-coordinate in 
    % "reference" data.
    if ~isempty(indexesXcoord)
        
        % Both X-coordinates match. Select all Y-coordinates of "dataToFilter".
        yCoordData = dataToFilter(indexesXcoord, 2);
        
        % Select Y-coordinates of "reference" data.
        yCoordRef = reference(reference(:, 1) == xCoordRef(iXCoord), 2);
        
        % Iterate over each Y-coordinate in "reference" data only in
        % X-coordinate previously selected
        for jYCoord = 1:size(yCoordRef)
            % First index where Y-coordinates of "dataToFilter" is equal to
            % Y-coordinate in "reference" data selected.
            firstIndex = find(yCoordData == yCoordRef(jYCoord), 1);
            
            % Save X-coordinate, Y-coordinate and Value.
            filteredData(kTrace, :) = dataToFilter(indexesXcoord(1) + firstIndex - 1, :);
            
            kTrace = kTrace + 1;
            
        end % End for

    end % End if
    
end % End for

% Remove zero rows
filteredData( ~any(filteredData, 2), : ) = [];