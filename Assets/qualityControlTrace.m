function qualityControlTrace(seismicTraces, params, coordinates, timeData, names, index)
% Make plots for quality control. Fractal Dimension.
%
% Input:
%  seismicTraces = Seismic traces selected.
%         params = Parameters: - Selected trace number.
%                              - Slope of the line log10(dividers) vs. log10(length)
%                              - Fractal Dimension.
%                              - R2.
%    coordinates = Coordinates of each trace selected.
%       timeData = Matrix with time windows for each calculation.
%          names = Cell array with seismic attributes names.
%          index = Trace index to plot.

% Select data
data = seismicTraces(:,:,index);

% Only plot if all elementes are non zero
if ~all(data(:)==0)

    figure(index)

    % 3D graph of two traces. Only use the first two seismic attributes 
    % used in Fractal Dimension calculation.

    traceNumber = num2str(params(index,1));
    xCoord = num2str(coordinates(index,1));
    yCoord = num2str(coordinates(index,2));
    
    titleFigure = ['Sesmic Trace #', traceNumber, '. Coordinates =  E:', xCoord, '  N:', yCoord];
    if size(seismicTraces,2)==1
        plot(data(:,1), timeData(:,index), 'b', 'LineWidth', 2);
        title(titleFigure);
        grid
        set(gca, 'YDir', 'reverse', 'LineWidth', 2, 'Color', 'w')
        axis('auto')
        xlabel(names{1})
        ylabel('Time')
    else
        plot3(data(:,1), data(:,2), timeData(:,index), 'b', 'LineWidth', 2);
        title(titleFigure);
        grid
        set(gca, 'ZDir', 'reverse', 'LineWidth', 1, 'Color', 'w', 'box', 'on')
        axis('vis3d') %normal
        xlabel(names{1});
        handleY= ylabel(names{2});
        zlabel('Time')

        % Change position of axis label
        handleY.Position=handleY.Position+[0.5 0 -1.5];
        
        % Rotate axis label
        set(get(gca,'xlabel'),'rotation', 20);
        set(get(gca,'ylabel'),'rotation',-35);
        
    end

end