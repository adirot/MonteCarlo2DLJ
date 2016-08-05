function plotParticles(coords,L,r,varargin)
            
% plots the particles with PBC

    p1 = inputParser();
    addOptional(p1, 'angles', []);
    addOptional(p1, 'dipoleStrength', []);
    parse(p1, varargin{:});
    Results1 = p1.Results;
    angles = Results1.angles;
    dipoleStrength = Results1.dipoleStrength;
    
    [~,N] = size(coords);
    figure; hold on;
    axis([-L/2 L/2 -L/2 L/2]);
    
    if isempty(dipoleStrength)
        dipoleStrength = ones(1,N);
    end
    
    % create color map:
    uniDipoleStrength = sort(unique(dipoleStrength));
    numOfColors = length(uniDipoleStrength);
    cmap = jet(numOfColors);
    
    for j = 1:N

        x = coords(1,j);
        y = coords(2,j);
        
        colorInd = (dipoleStrength(j) == uniDipoleStrength);
        color = cmap(colorInd,:);
        
        if ~isempty(angles)
            plotCircle(x,y,r,angles(1,j),'color',color);
        else
            plotCircle(x,y,r,'color',color);
        end

        %PBC
        if  x > (L/2 - r)
            if ~isempty(angles)
                plotCircle(x-L,y,r,angles(1,j),'color',color);
            else
                plotCircle(x-L,y,r,'color',color);
            end

        else if x < (-L/2 + r)
            if ~isempty(angles)
                plotCircle(x+L,y,r,angles(1,j),'color',color);
            else
                plotCircle(x+L,y,r,'color',color);
            end

            end
        end

        if  y > (L/2 - r)

            
            if ~isempty(angles)
                plotCircle(x,y-L,r,angles(1,j),'color',color);
            else
                plotCircle(x,y-L,r,'color',color);
            end

        else if y < (-L/2 + r)

            if ~isempty(angles)
                plotCircle(x,y+L,r,angles(1,j),'color',color);
            else
                plotCircle(x,y+L,r,'color',color);
            end

            end
        end

        % plot corner particles

        if  x > (L/2 - r) && y > (L/2 - r)
            if ~isempty(angles)
                plotCircle(x-L,y-L,r,angles(1,j),'color',color);
            else
                plotCircle(x-L,y-L,r,'color',color);
            end

        end

        if  x > (L/2 - r) && y < (-L/2 + r)
            if ~isempty(angles)
                plotCircle(x-L,y+L,r,angles(1,j),'color',color);
            else
                plotCircle(x-L,y+L,r,'color',color);
            end

        end

        if  x < (-L/2 + r) && y > (L/2 - r)
            if ~isempty(angles)
                plotCircle(x+L,y-L,r,angles(1,j),'color',color);
            else
                plotCircle(x+L,y-L,r,'color',color);
            end

        end

        if  x < (-L/2 + r) && y < (-L/2 + r)
            if ~isempty(angles)
                plotCircle(x+L,y+L,r,angles(1,j),'color',color);
            else
                plotCircle(x+L,y+L,r,'color',color);
            end

        end

    end

% Add legend
for i = 1:length(uniDipoleStrength)
    h(i) = plot(0,0,'color',cmap(i,:), 'visible', 'off');
    leg{1,i} = num2str(uniDipoleStrength(i));
end
l = legend(h, leg);
t = get(l,'title');
set(t,'string','Dipole strength');

    function plotCircle(x,y,r,varargin)
        p = inputParser();
        addOptional(p, 'angle', []);
        addOptional(p, 'color', 'b');
        parse(p, varargin{:});
        Results = p.Results;
        angle = Results.angle;
        col = Results.color;
        
        % plot centers
        plot(x,y,'+r');
        
        % plot boundarys
        ang=0:0.01:2*pi; 
        xp=r*cos(ang);
        yp=r*sin(ang);
        plot(x+xp,y+yp,'color',col);
        
        % plot oriantation
        if ~isempty(angle)
            plot([x-r*cos(angle) x+r*cos(angle)],...
                [y-r*sin(angle) y+r*sin(angle)],'r');
        end
    end
end
