function U = pairU(dist,rCutoff,m,varargin)

% calculates the reduced energy according to the pair
% potantial, only pair closer than rCutoff are regarded.

% uses the potantial:
% 4*(((1/r)^12)-((1/r).^m))

% input: dist is a row vector of all pair distances
% output: U is the total energy 

% the potantial can have angle dependence, given by an anonymus function
% of two angles. for a pair of two cells, A and B, the first angle is the
% oriantation of cell B with respect to the oriantation of A,
% and the second angle is the angle of the vector connecting the two cells,
% with respect to the oriantation of A. 
% for example, an angle dependence might be:
% f(a,b) = @(a,b) cos(a)*sin(b)
% meaning that the distance dependent potantial is multiplied by this
% factor.
% so we will calculate the pair potantial with:
% U = pairU(dist,rCutoff,m,'angleDependence',f,'relativeCellAngles',relativeCellAngles)
% where relativeCellAngles is a vactor of 2 by N, noteing the two angles
% needed as described.

    p = inputParser();
    addOptional(p, 'angleDependence', []);
    addOptional(p, 'relativeCellAngles', []);
    parse(p, varargin{:});
    Results = p.Results;
    angleDependence = Results.angleDependence;
    relativeCellAngles = Results.relativeCellAngles;
    
    cutoffInd = dist < rCutoff;
    dist_lt_rCutoff = dist(cutoffInd);
    % u is the energies of each pair
    u = 4*(((1./dist_lt_rCutoff).^12)-((1./dist_lt_rCutoff).^m)); 
    
    if ~isEmpty(angleDependence)
        angleFactor = angleDependence(...
            relativeCellAngles(1,:),relativeCellAngles(2,:));
        u = u.*angleFactor;
    end
    
    U = sum(u);
end