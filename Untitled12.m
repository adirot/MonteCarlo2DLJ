numOfSquares = (1:30).^2;
L = sqrt(625/0.005);

for i = 1:length(numOfSquares)
    [~, ~, ns] =...
        rhoDistribution2(coords{1,1,2},L,numOfSquares(i),10);
    numOfCellsInSquares{1,i} = ns;
    rho{1,i} = min(numOfCellsInSquares{i}):max(numOfCellsInSquares{i});
    [~, PL{1,i}, ~] =...
        rhoDistribution2(coords{1,1,2},L,numOfSquares(i),length(rho{1,i}));
    
    if i == 1 
        maxLenRho = length(rho{1,i});
    else
        if maxLenRho < length(rho{1,i})
            maxLenRho = length(rho{1,i});
        end
    end
end    

x = zeros(length(numOfSquares),maxLenRho);
y = zeros(length(numOfSquares),maxLenRho);

for i = 1:length(numOfSquares)
    x(i,1:length(rho{1,i})) = rho{1,i};
    y(i,1:length(rho{1,i})) = PL{1,i}.*rho{1,i}/625;
    
    length2plot(i) = length(rho{1,i});
    legend{1,i} = ['number of squares: ' num2str(numOfSquares(i))];

end

colorPlot(x,y,'length2plot',length2plot,'addLegend',legend);