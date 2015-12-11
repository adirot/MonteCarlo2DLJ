function dist = distPBC(x,y,allPossition,L)
            % calculate the PBC distances of the point (x,y) from all
            % particle possitions
            
                distx = abs(x - allPossition(1,:));
                disty = abs(y - allPossition(2,:));
                bigDistx = distx > (L/2);
                bigDisty = disty > (L/2);
                distx = distx - L.*bigDistx;
                disty = disty - L.*bigDisty;
                dist = sqrt(distx.^2 + disty.^2);
end