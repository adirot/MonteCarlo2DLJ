figure;
plot(coords{1,1,2}(1,:),coords{1,1,2}(2,:),'+');
boxSide = L;
numOfSquares = 900;
xlim([-boxSide/2 boxSide/2]);
ylim([-boxSide/2 boxSide/2]);
hold on;
squareSide = (boxSide/sqrt(numOfSquares));
m = sqrt(numOfSquares);
for i = (-boxSide/2):squareSide:(boxSide/2 - squareSide)
    plot([i i],[-boxSide/2 boxSide/2],'r');
    plot([-boxSide/2 boxSide/2],[i i],'r');
end
