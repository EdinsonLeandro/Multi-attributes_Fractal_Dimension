function qualityControlDividerLength(seismicTraces, params, dividersLength, index)
% Make plots for quality control. Fractal Dimension.
%
% Input:
%  seismicTraces = Seismic traces selected.
%         params = Parameters: - Selected trace number.
%                              - Slope of the line log10(dividers) vs. log10(length)
%                              - Fractal Dimension.
%                              - R2.
% dividersLength = Dividers and Total Length. 3rd Dimension is each trace.
%          index = Trace index to plot.

% Select data
data = seismicTraces(:,:,index);

% Only plot if all elementes are non zero
if ~all(data(:)==0)

    figure(index)

    % Plot: Fractal Dimension.
    slope = num2str(params(index,2));
    traceNumber = num2str(params(index,1));

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
    % xlim = get(gca,'XLim');
    % ylim = get(gca,'YLim');
    %  text(mean(xlim), 1.4*ylim(2),['R^2 = ', num2str(params(index,4))],...
    %      'FontWeight','bold','HorizontalAlignment','center');
    ylabel('Log10 (Total Length)'), xlabel('Log10 (Divider)');
    title(['Trace #', traceNumber, '.  Slope=', slope, '.  R^2=', ...
        num2str( round(params(index,4), 2) )]);

end