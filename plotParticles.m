function plotParticles(coords,L,r,varargin)
            
% plots the particles with PBC

    p = inputParser();
    addOptional(p, 'angles', []);
    parse(p, varargin{:});
    Results = p.Results;
    angles = Results.angles;

    [~,N] = size(coords);
    figure; hold on;
    axis([-L/2 L/2 -L/2 L/2]);

    for j = 1:N

        x = coords(1,j);
        y = coords(2,j);
        plotCircle(x,y,r,angles(1,j));

        %PBC
        if  x > (L/2 - r)

            plotCircle(x-L,y,r,angles(1,j));

        else if x < (-L/2 + r)

                plotCircle(x+L,y,r,angles(1,j));

            end
        end

        if  y > (L/2 - r)

            plotCircle(x,y-L,r,angles(1,j));

        else if y < (-L/2 + r)

            plotCircle(x,y+L,r,angles(1,j));

            end
        end

        % plot corner particles

        if  x > (L/2 - r) && y > (L/2 - r)

            plotCircle(x-L,y-L,r,angles(1,j));

        end

        if  x > (L/2 - r) && y < (-L/2 + r)

            plotCircle(x-L,y+L,r,angles(1,j));

        end

        if  x < (-L/2 + r) && y > (L/2 - r)

            plotCircle(x+L,y-L,r,angles(1,j));

        end

        if  x < (-L/2 + r) && y < (-L/2 + r)

            plotCircle(x+L,y+L,r,angles(1,j));

        end

    end



    function plotCircle(x,y,r,varargin)
        p = inputParser();
        addOptional(p, 'angle', []);
        parse(p, varargin{:});
        Results = p.Results;
        angle = Results.angle;
        
        % plot centers
        plot(x,y,'+r');
        
        % plot boundarys
        ang=0:0.01:2*pi; 
        xp=r*cos(ang);
        yp=r*sin(ang);
        plot(x+xp,y+yp);
        
        % plot oriantation
        if ~isempty(angle)
            plot([x-r*cos(angle) x+r*cos(angle)],...
                [y-r*sin(angle) y+r*sin(angle)],'r');
        end
    end
end
