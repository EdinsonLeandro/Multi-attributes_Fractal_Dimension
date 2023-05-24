function qualityControlHist(data)
% R^2 histogram in Fractal Dimension Calculation.

figure('Name', 'R^2 in Fractal Dimension Calculation')

% Remove zero rows
data( ~any(data, 2), : ) = [];

% Plot
hist(data, (0.1:0.1:0.9));
ylabel('Frequency'), xlabel('R^2');
title('R^2 in Fractal Dimension Calculation');