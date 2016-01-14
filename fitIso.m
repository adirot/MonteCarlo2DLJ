function fit = fitIso(isotherms,fitprop,varargin)
        
    fit = [];

     p = inputParser();
     addOptional(p, 'figureHandle', []);
     addOptional(p, 'my_eps', 0.1); % resulution for the virial
     addOptional(p, 'colorVir', false);
     parse(p, varargin{:});
     Results = p.Results;
     figureHandle = Results.figureHandle;
     my_eps = Results.my_eps;
     colorVir = Results.colorVir;
    
     for i = 1:length(isotherms)
         switch fitprop
             case 'linearFit'
                    [fit{i,1}, fit{i,2}] = ...
                        polyfit(isotherms(i).rho,isotherms(i).pressure,1);
             case 'plotVirialExp'
                        
                 [P,rho,~,good_ind,cantusevirial_ind,bad_ind] = ...
                    real_pressure2D(isotherms(i).T,...
                            isotherms(i).rho,my_eps,'');
                fit(i).P = P;
                fit(i).rho = rho;
                fit(i).good_ind = good_ind;
                fit(i).cantusevirial_ind = cantusevirial_ind;
                fit(i).bad_ind = bad_ind;
                
         end
     end
     
     if ~isempty(figureHandle)
         figure(figureHandle);
         hold on;
         x = linspace(min(isotherms(1).rho),...
                     max(isotherms(1).rho),100);
                 
         switch fitprop
             
             case 'linearFit'
                 
                 for i = 1:length(isotherms)
                    p = fit{i,1};
                    y = p(1)*x + p(2);
                    plot(x,y);
                    text(x(70),y(70),['T = ' num2str(p(1))...
                        ' (' num2str(fit{i,2}.R(1,1)) ','...
                        num2str(fit{i,2}.R(1,2)) ')']);
                    fit(i).P = y; 
                 end
                 
             case 'plotLin'
                 dim = [.2 .5 .3 .3];
                 str = 'P = T\cdot\rho'; 
                 annotation('textbox',dim,'String',str,...
                     'FitBoxToText','on','Color','red');
                    
                 for i = 1:length(isotherms)
                    
                    plot(isotherms(i).rho,isotherms(i).T*isotherms(i).rho,'r');
                    fit(i).P = isotherms(i).T*isotherms(i).rho;
                    
                 end 
                 
             case 'plotVirialExp'
                 dim = [.2 .5 .3 .3];
                 str = 'P = T\cdot\rho (1 - \rho\cdot B_2)'; 
                 annotation('textbox',dim,'String',str,...
                     'FitBoxToText','on');
                 
                 if colorVir
                     for i = 1:length(isotherms)

                        plot(fit(i).rho(good_ind),fit(i).P(good_ind),...
                            'color',[0 0.5 0],'marker','o','line','none'); %green

                     end

                     for i = 1:length(isotherms)

                        plot(fit(i).rho(cantusevirial_ind),...
                            fit(i).P(cantusevirial_ind),...
                            'color','y','marker','o','line','none'); 

                     end


                     for i = 1:length(isotherms)

                        plot(fit(i).rho(bad_ind),fit(i).P(bad_ind),...
                            'color','r','marker','o','line','none'); 

                     end
                 else
                     for i = 1:length(isotherms)

                        plot(fit(i).rho,fit(i).P,...
                            'color','k','marker','o','line','none'); 

                     end
                     
                 end
                 
         end
     end
            
            
    
end