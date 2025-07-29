function StatusBar(currentIdx, totalIdx, label, skip)
    %	 ____________________________________________
    %	|                                            |
    %	|     S  T  A  T  U  S       B  A  R         |
    %	|____________________________________________|
    %	| |~~|  |~~|  |~~|  |~~|  |~~|  |~~|  |~~|    |
    %	| |__|  |__|  |__|  |__|  |__|  |__|  |__|    |
    %	|                                            |
    %	|   Welcome to the Status Bar!               |
    %	|   - Now serving real-time updates          |
    %	|   - Progress on tap!                       |
    %	|____________________________________________|
    %	 
    %	     |   |   |   |   |   |   |   |   |   |  
    %	     |___|___|___|___|___|___|___|___|___|   
    %	Displays an in-place updating progress bar in the MATLAB Command Window.
    %	Only updates every 'skip' iterations (optional).
    %
    %	Inputs:
    %	  currentIdx  - Current loop index (e.g., i)
    %	  totalIdx    - Total number of iterations (e.g., N)
    %	  label       - (Optional) String label, e.g., 'File', 'Image'
    %	  skip        - (Optional) Update frequency (e.g., skip=10 updates every 10th iteration)
    %	Author: Felix Maurer
    %	License: MIT

    
    %	Handle optional inputs
    if nargin < 3
        label = '';
    end
    if nargin < 4
        skip = 1;
    end

    % Only update on multiples of 'skip' or at the final iteration
    if mod(currentIdx, skip) ~= 0 && currentIdx ~= totalIdx
        return;
    end

    % Persistent variable to track previous output length
    persistent lastLength;
    if isempty(lastLength) || currentIdx == 1
        lastLength = 0;
    end

    % Compute progress bar string
    percentDone = currentIdx / totalIdx;
    barLength = 30;
    numBars = round(percentDone * barLength);
    barStr = [repmat('#', 1, numBars), repmat('-', 1, barLength - numBars)];

    if ~isempty(label)
        suffix = sprintf('%s %d of %d', label, currentIdx, totalIdx);
    else
        suffix = sprintf('%d of %d', currentIdx, totalIdx);
    end

    statusStr = sprintf('[%s] %3.0f%% — %s', ...
        barStr, percentDone * 100, suffix);

    % Erase previous output
    fprintf(repmat('\b', 1, lastLength));
    fprintf('%s', statusStr);

    % Store new output length
    lastLength = length(statusStr);

    drawnow;

    % Clean up at end
    if currentIdx == totalIdx
        fprintf('\n');
        lastLength = 0;
    end
end
