function plotParticles(coords,L,r)
            
            % plots the particles with PBC
            
                [~,N] = size(coords);
                figure; hold on;
                axis([-L/2 L/2 -L/2 L/2]);
                
                for j = 1:N
                    
                    x = coords(1,j);
                    y = coords(2,j);
                    plotCircle(x,y,r);
                    
                    %PBC
                    if  x > (L/2 - r)
                        
                        plotCircle(x-L,y,r);
                        
                    else if x < (-L/2 + r)
                            
                            plotCircle(x+L,y,r);
                          
                        end
                    end
                                        
                    if  y > (L/2 - r)
                        
                        plotCircle(x,y-L,r);
                        
                    else if y < (-L/2 + r)
                           
                        plotCircle(x,y+L,r);
                        
                        end
                    end
                    
                    % plot corner particles

                    if  x > (L/2 - r) && y > (L/2 - r)
                        
                        plotCircle(x-L,y-L,r);
                        
                    end
                    
                    if  x > (L/2 - r) && y < (-L/2 + r)
                        
                        plotCircle(x-L,y+L,r);
                        
                    end
                    
                    if  x < (-L/2 + r) && y > (L/2 - r)
                        
                        plotCircle(x+L,y-L,r);
                        
                    end
                    
                    if  x < (-L/2 + r) && y < (-L/2 + r)
                        
                        plotCircle(x+L,y+L,r);
                        
                    end

                end
                
                
                                
            function plotCircle(x,y,r)
                        % plot centers
                            plot(x,y,'+r');
                        % plot boundarys
                            ang=0:0.01:2*pi; 
                            xp=r*cos(ang);
                            yp=r*sin(ang);
                            plot(x+xp,y+yp);
            end
        end
