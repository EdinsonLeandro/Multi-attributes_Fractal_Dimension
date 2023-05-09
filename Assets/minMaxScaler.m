function X=minMaxScaler(data, values)
% X = minMaxScaler(data, values) performs Min Max Scaler normalization:
% Xscaled = (X - Xmin) / (Xmax - Xmin)
%
% Input parameters:
% data: Array with the data to be normalized.
% values: Matrix with minimum and maximun values.
%
% Output Parameters:
% X: Array with each normalized data

% Assumption: The number of columns of "values.'" and "data" are equal.
% Transpose matrix with minimun and maximun values
values = values.';

% Initialize "X" matrix with normalized data
X = zeros(size(data));

for iColumn = 1:size(data,2)
    minValue = values(1, iColumn);
    maxValue = values(2, iColumn);

    X(:, iColumn) = ( data(:, iColumn)-minValue) / (maxValue-minValue);
end