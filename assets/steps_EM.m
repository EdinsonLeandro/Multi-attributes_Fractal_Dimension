function nSteps = steps_EM(X, divider)
% The original function: n = steps3(X, d)
% calculates the number of steps to walk along the set of points X
% using a divider length d. 
% Written by Gerry Middleton, McMaster Univ., June 1996
% EMail: middleto@mcmaster.ca
% calls nextpt3.m To increase speed at the risk of accuracy
% change this to nextpt.
%
% This is a modified version to calculate the steps number
% needed to approximate the trace X with N-Dimension with "divider".

nSteps = 0;
index = 2;

% Xn = From this point the distance is measured to the next sample.
% The algotithm starts with the first point
Xn=X(1,:);

[nSamples, ~] = size(X);
while index <= nSamples
    % It will compare with the given divisor.

    % Euclidean distance to the next sample.
    diff = Xn - X(index,:);
    distanceNext = sqrt(sum(diff.^2,2));
    
    % Euclidean distance to previous sample.
    diff = Xn - X(index-1, :);
    distancePrevious = sqrt(sum(diff.^2,2));
    
    % di = distanceNext
    % dim = distancePrevious

    if distanceNext >divider
        if distancePrevious == 0
            ratio = distanceNext/divider;
            nSteps = nSteps+ratio;
            
            % Select next point
            Xn= X(index,:);
            index = index + 1;
        else
            % distancePrevious ~= 0
            dk = distanceNext - divider;
            dj = divider - distancePrevious;
            if dk>dj
                Xn= X(index-1,:);
                nSteps= nSteps+1;
            else
                Xn = X(index,:);
                nSteps = nSteps + 1;
                index = index + 1;
            end
            
            if dk == dj
                ratio = distancePrevious/divider;
                nSteps = nSteps+ratio;
                Xn = X(index-1,:);
            end
        end
    
    elseif distanceNext < divider
        if index == nSamples
            ratio = distanceNext/divider;
            nSteps = nSteps+ratio;
            break
        end
        %         
        % For distances less than "divider", the fraction "ratio" is calculated
        %         
        index=index+1;
    
    else
        Xn = X(index,:);
        nSteps = nSteps+1;
        index = index+1;
    end
    
end
