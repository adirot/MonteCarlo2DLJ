function hist = my_hist(data,bins)
% makes a histogram from the data. x is the x axis of the histogram.
% x must be a linearly spaced vector.

% to plot the result:
 
% increment = bins(2)-bins(1);
% bar(bins+increment/2,hist,1);

    % find the bin index for every data point
    binIndex4data = ceil((data-bins(1))/(bins(2) - bins(1)));
    binIndex = 1:length(bins);
    
    % equals(i,j) = 1 means that binIndex4data(i) equals binIndex(j).
    equals = bsxfun(@eq,binIndex4data',binIndex);
    
    % counts how many data point are in each bin (sums every cloumn of
    % equals matrix) 
    hist = sum(equals);
    
end