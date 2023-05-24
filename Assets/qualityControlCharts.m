function qualityControlCharts(seismicTraces, params, dividersLength, coordinates, timeData, names, index)
% Make plots for quality control. Fractal Dimension.
%
% Input:
%  seismicTraces = Seismic traces selected.
%         params = Parameters: - Selected trace number.
%                              - Slope of the line log10(dividers) vs. log10(length)
%                              - Fractal Dimension.
%                              - R2.
% dividersLength = Dividers and Total Length. 3rd Dimension is each trace.
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
        subplot(1,2,1);
        plot(data(:,1), timeData(:,index), 'b', 'LineWidth', 2);
        title(titleFigure);
        grid
        set(gca, 'YDir', 'reverse', 'LineWidth', 2, 'Color', 'w')
        axis('auto')
        xlabel(names{1})
        ylabel('Time')
    else
        subplot(1,2,1);
        plot3(data(:,1), data(:,2), timeData(:,index), 'b', 'LineWidth', 2);
        title(titleFigure);
        grid
        set(gca, 'ZDir', 'reverse', 'LineWidth', 1, 'Color', 'w', 'box', 'on')
        axis('vis3d') %normal
        % handleX= xlabel(names{1});
        handleY= ylabel(names{2});
        zlabel('Time')

        % Change position of axis label
        % handleX.Position=handleX.Position+[-10.0 -12.0 -2.5];
        handleY.Position=handleY.Position+[0.5 0 -3.0];
        
        % Rotate axis label
        set(get(gca,'xlabel'),'rotation', 20);
        set(get(gca,'ylabel'),'rotation',-35);
        
    end

    % Plot: Fractal Dimension.
    slope = num2str(params(index,2));

    % Select "dividers" and "length"
    dividers = dividersLength(1,:,index);
    length = dividersLength(2,:,index);
    
    % Remove zeros rows (if exist) within the matrix
    % (which are generated due to variability of matrix dimensions).
    dividers( ~any(dividers, 2), : ) = [];
    length( ~any(length, 2), : ) = [];
  
    % Logarithm base 10, for divider and total length. 
    logDividers = log10(dividers);
    logLength = log10(length);

    % Linear regression
    coefficients = polyfit(logDividers,logLength,1);
    yValues = polyval(coefficients,logDividers);

    hold on
    subplot(1,2,2);
    plot(logDividers, logLength, 'ks', logDividers, yValues, 'g',...
        'LineWidth', 2, 'MarkerFaceColor', 'b');
    axis('equal');
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    text(mean(xlim), 1.4*ylim(2),['R^2 = ', num2str(params(index,4))],...
        'FontWeight','bold','HorizontalAlignment','center');
    ylabel('Log10 (Total Length)'), xlabel('Log10 (Divider)');
    title(['Seismic Trace #' traceNumber '.     Slope = ' slope]);

end