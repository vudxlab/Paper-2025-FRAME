function M = TriDiagonalMatrix(n)
    % Check if the input is less than 3
    if n < 3
        M = NaN;
        return;
    end
    
    % Pre-allocate the matrix M with zeros
    % Pre-allocate the matrix M with zeros
    M = zeros(n);
    
    % Fill the middle rows according to the pattern
    for i = 2:n-1
        M(i, i-1) = i ; % Fill the sub-diagonal with i-1
        M(i, i)   = i * 2; % Fill the main diagonal with i*2
        M(i, i+1) = i * 3; % Fill the super-diagonal with i*3
    end
end
