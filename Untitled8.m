for i = (-L/2):squareSide:(L/2 - squareSide)
    indi = indi + 1;
    indj = 0;
    for j = (-L/2):squareSide:(L/2 - squareSide)
        indj = indj + 1;
        densities(indi,indj) = countParticlesInSquare(i,j,coords(1,:),coords(2,:));
    end
end


   function N = countParticlesInSquare(i,j,x,y)
    % Count the number of particles in square i,j
        N = sum(and(and(y > j,y < j+1), and(x > i,x < i+1)));
    end