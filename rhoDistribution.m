function [rho, PL] = rhoDistribution(coords,L,numOfSquares,numOfBins)
%% Find a histogram of the densities in a montecarlo step

% Given coordinates of N particles in 2D box, we devide the box to numOfSquares
% squares and build a histogram from the densities in each square.

% Input:
% ~~~~~
% coords       : Coordinates in 2D, 2 by N matrix. the center of the box
%                must be at (0,0);
% L            : Size of box side.
% numOfSquares : Number of squares we devide the board into. the root of
%                this number must be real
% numOfBins    : Number of bins in the histogram

% Output:
% ~~~~~
% rho : X axis of the histogram, densities.
% PL  : Y axis of the histogram, the number of densities found for each
%       bin, normalized by rho0.

% Usage example:
% ~~~~~~~~~~~~~
% % create N random coordinates in a boxSide by boxSide box.
% N = 100;
% boxSide = 1;
% coords(1,:) = boxSide*rand(1,N) - boxSide/2;
% coords(2,:) = boxSide*rand(1,N) - boxSide/2;

% % get the density distribution:
% numOfSquares = 100;
% [rho, PL] = rhoDistribution(coords,boxSide,numOfSquares,10)

% % plot results

% % plot particles and squares 
% plot(coords(1,:),coords(2,:),'+');
% xlim([-boxSide/2 boxSide/2]);
% ylim([-boxSide/2 boxSide/2]);
% hold on;
% squareSide = sqrt(squareArea);
% m = sqrt(numOfSquares);
% for i = (-boxSide/2):squareSide:(boxSide/2 - squareSide)
%     plot([i i],[-boxSide/2 boxSide/2],'r');
%     plot([-boxSide/2 boxSide/2],[i i],'r');
% end

% % plot histogram
% figure; plot(rho,PL);

% Get the density in each square 

m = sqrt(numOfSquares);
if isinteger(m)
    error('The square root of numOfSquares must be an intiger');
end

squareArea = L^2/numOfSquares;
squareSide = sqrt(squareArea);
densities = zeros(m-1,m-1);

for i = (-boxSide/2):squareSide:(boxSide/2 - squareSide)
    for j = (-boxSide/2):squareSide:(boxSide/2 - squareSide)
        densities(i,j) = countParticlesInSquare(i,j,coords(1,:),coords(2,:));
    end
end

densities = densities/squareArea;

% Bin the results to a histogram
hist(densities,numOfBins);



%% Functions used in this code:

    function N = countParticlesInSquare(i,j,x,y)
    % Count the number of particles in square i,j
        N = sum(and(and(y > j,y < j+1), and(x > i,x < i+1)));
    end
end