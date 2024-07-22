function indices = getCircleIndices(centerIndex, radius, grid)
    % Validate inputs
    if numel(centerIndex) ~= 2
        error('centerIndex must be a 2-element vector for a 2D grid.');
    end
    if radius <= 0
        error('Radius must be positive.');
    end
    
    % Get grid dimensions
    [numRows, numCols] = size(grid);
    
    % Create coordinate grids
    [rowGrid, colGrid] = ndgrid(1:numRows, 1:numCols);
    
    % Calculate squared distances from the central point
    distSquared = (rowGrid - centerIndex(1)).^2 + (colGrid - centerIndex(2)).^2;
    
    % Find indices within the circle
    withinCircle = distSquared <= radius^2;
    
    % Extract row and column indices within the circle
    [rowIndices, colIndices] = find(withinCircle);
    
    % Ensure indices are within bounds
    validRows = rowIndices(rowIndices >= 1 & rowIndices <= numRows);
    validCols = colIndices(colIndices >= 1 & colIndices <= numCols);
    
    % To maintain matching rows and columns:
    [validRows, validCols] = deal(validRows(ismember(rowIndices, validRows)), ...
                                   validCols(ismember(colIndices, validCols)));
    
    % Combine row and column indices into a cell array
    cellIndices = {validRows, validCols};

    %
    indices = sub2ind(size(grid), cellIndices{1}, cellIndices{2});
end