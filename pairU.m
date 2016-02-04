function U = pairU(dist,rCutoff,m)

% calculates the reduced energy according to the pair
% potantial, only pair closer than rCutoff are regarded.

% input: dist is a row vector of all pair distances
% output: U is the total energy 

    dist_lt_rCutoff = dist(dist < rCutoff);
    % u is the energies of each pair
    u = 4*(((1./dist_lt_rCutoff).^12)-((1./dist_lt_rCutoff).^(-m))); 
    U = sum(u);
end