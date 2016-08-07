function U = pairU(dist,rCutoff,m,varargin)

% calculates the reduced energy according to the pair
% potantial, only pair closer than rCutoff are regarded.

% uses the potantial:
% (((1/r)^12)-((1/r).^m))

% input: dist is a row vector of all pair distances
% output: U is the total energy 

% the potantial can have angle dependence, given by an anonymus function
% of two angles. for a pair of two cells, A and B, the first angle is the
% oriantation of cell B with respect to the oriantation of A,
% and the second angle is the angle of the vector connecting the two cells,
% with respect to the oriantation of A. 
% for example, an angle dependence might be:
% f(a,b) = @(a,b) cos(a)*cos(b)
% meaning that the distance dependent potantial is multiplied by this
% factor. 
% so we will calculate the pair potantial with:
% U = pairU(dist,rCutoff,m,'angleDependence',f,'relativeCellAngles',relativeCellAngles)
% where relativeCellAngles is a vactor of 2 by N, noteing the two angles
% needed as described.

    p = inputParser();
    addOptional(p, 'angleDependence', []);
    addOptional(p, 'relativeCellAngles', []);
    addOptional(p, 'numOfrelativeCellAngles', 2);
    addOptional(p, 'ufunc', @(r) (((1./r).^12)-((1./r).^m)));
    addOptional(p, 'dipolePairStrength', []);
    parse(p, varargin{:});
    Results = p.Results;
    angleDependence = Results.angleDependence;
    relativeCellAngles = Results.relativeCellAngles;
    numOfrelativeCellAngles = Results.numOfrelativeCellAngles;
    ufunc = Results.ufunc;
    dipolePairStrength = Results.dipolePairStrength;
    
    cutoffInd = dist < rCutoff;
    dist_lt_rCutoff = dist(cutoffInd);
    % u is the energies of each pair
    u = ufunc(dist_lt_rCutoff); 
    if isempty(dipolePairStrength)
        dipolePairStrength = ones(1,length(dist_lt_rCutoff));
    else
        dipolePairStrength = dipolePairStrength(cutoffInd);
    end
    
    if isempty(u)
        u = 0;
    else
        if ~isempty(angleDependence)
            if numOfrelativeCellAngles == 2
                angleFactor = angleDependence(...
                    relativeCellAngles(cutoffInd,1),relativeCellAngles(cutoffInd,2));
            else
                if numOfrelativeCellAngles == 3
                    angleFactor = angleDependence(...
                    relativeCellAngles{1,1}(cutoffInd),...
                    relativeCellAngles{1,2}(cutoffInd),...
                    relativeCellAngles{1,3}(cutoffInd));
                else
                    error('numOfrelativeCellAngles must be 2 or 3');
                end
            end
            

            u = (u.*angleFactor).*dipolePairStrength;
        end
    end
    
    U = sum(u);
end
