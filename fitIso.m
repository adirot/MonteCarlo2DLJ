function fit = fitIso(isotherms,fitprop,varargin)
        
    fit = [];

     p = inputParser();
     addOptional(p, 'figureHandle', []);
     parse(p, varargin{:});
     Results = p.Results;
     figureHandle = Results.figureHandle;
    
     for i = 1:length(isotherms)
         switch fitprop
             case 'linearFit'
                    [fit{i,1}, fit{i,2}] = ...
                        polyfit(isotherms(i).rho,isotherms(i).pressure,1);
             case 'plotVirialExp'
                 [P,rho,~,~,~,~] = ...
                    real_pressure2D(isotherms(i).T,isotherms(i).rho,0.1,'');
                fit(i).P = P;
                fit(i).rho = rho;
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
                 end
                 
             case 'plotLin'
                 dim = [.2 .5 .3 .3];
                 str = 'P = T\cdot\rho'; 
                 annotation('textbox',dim,'String',str,...
                     'FitBoxToText','on','Color','red');
                    
                 for i = 1:length(isotherms)
                    
                    plot(x,isotherms(i).T*x,'r');
                    
                 end 
                 
             case 'plotVirialExp'
                 dim = [.2 .5 .3 .3];
                 str = 'P = T\cdot\rho (1 - \rho\cdot B_2)'; 
                 annotation('textbox',dim,'String',str,...
                     'FitBoxToText','on');
                 
                 for i = 1:length(isotherms)
                    
                    plot(fit(i).rho,fit(i).P,'ok');
                    
                 end 
                 
         end
     end
            
            
    
end