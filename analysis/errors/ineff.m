    function tauProp = ineff(prop,n,firstSteps2ignore)
           % for error calculations. see computer simulation of liquids
           % page 192
           
           
           % calculate <A>
           meanProp = mean(prop);
           
           prop = prop(firstSteps2ignore:length(prop));
           
           nt = (length(prop) - firstSteps2ignore);
           
           for j = 1:length(n)
                % calculate <A>b
                ind = 1;
                for i = 1:n(j):(nt-n(j)+1)
                   
                    Pmeanb(ind) = mean(prop(i:(i+n(j)-1)));
                    ind = ind + 1;
                end
           
                % calculate sigma^2(<A>b) for all the different tau values
                
                nb(j) = length(Pmeanb);
                varMeanProp(j) = mean((Pmeanb - mean(Pmeanb)).^2);
                Pmeanb = [];
                
           end
            
           %calculate sigma^2(A)
           varProp = mean((prop(:) - meanProp).^2);
           
           % calculate tau
           tauProp = n.*varMeanProp/varProp;
            
           
               figure;
               plot(sqrt(n),tauProp);
               hold on;
               title('s for prop');
               xlabel('\sqrt{n}');
               ylabel('s');
           
           
       end
   