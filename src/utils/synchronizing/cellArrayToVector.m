function [outputVector] = cellArrayToVector(cellArray)
% This function literally just takes all of the values in all the cells of
% a cell array and outputs one column vector with them.
%
% Adam Smoulder, 7/9/20

% Flatten cell array
cellArray = cellArray(:);

% Loop over cells and get all values
outputVector = [];
for i = 1:length(cellArray)
    outputVector = [outputVector ; cellArray{i}(:)];
end; clear i

end

