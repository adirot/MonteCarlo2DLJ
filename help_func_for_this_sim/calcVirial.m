function virial = calcVirial(dists,rho,n,m,N,rCutoff)

%% calculate the reduced virial of a system of 2D particles with pair potantial
%       using the distances between particles 'allDist',
%       calculate the virial.

%        ==================================
%       | virial = - (1/2)sum(r_ij (dU/dr)) |
%        ==================================

%       in reduced units, with a 4epsilon[(1/r)^n - (1/r)^m] potantial
%       this is equivivalent to:
%    
%        ===========================================================
%       | virial' = - (2*rho'/N)*sum(m(1/r'_ij)^m - n(1/r'_ij)^n) |
%        ===========================================================

%       where r_ij is the distance between particles i,j. the sum is on all
%       particle pairs, each pair is counted once.
%       for more information about the first equation: http://www2.msm.ctw.utwente.nl/sluding/TEACHING/APiE_Script_v2011.pdf
%       eq. (3.15)       
%       for more information about the second equation, read 'a note about
%       reduced units' in the git repository

%       input: 
%       allDist - a matrix containing all pair distances for each step in a
%                 monte carlo simulation, created with 'MonteCarlo2DLJ'
%       rho - density in reduced units
%       n - the pair potential power law in the distance 
%           (pair potential in reduced units is u = 4*[(1/r)^n - (1/r)^m] 
%           where r is the pair distance)
%       m - the 'm' constant in the pair potantial
%       T - reduced Temperature 

       dists = dists(and(dists < rCutoff,dists > 0));
       
       % V = - (2*rho'/N)*sum(m(1/r'_ij)^m - n(1/r'_ij)^n)
       vir = - (2*rho/N)*(sum(sum(m*dists.^(-m) - n*dists.^(-n))));
       % make vir a row vector
       virial(1,:) = vir(1,1,:);

end