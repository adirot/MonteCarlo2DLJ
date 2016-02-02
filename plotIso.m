function [isotherms,fit,canGetUfromgRind,hPvsRho,hPvsV,U,T] = plotIso(varargin)
    
    % check input, build isotherm object from data files if necessary

        p = inputParser();
        addOptional(p, 'N', []); 
        addOptional(p, 'Visible', 'on');
        addOptional(p, 'saveFig', true);
        addOptional(p, 'fitprop', []);
        addOptional(p, 'residuals', []);
        addOptional(p, 'fileNameEnd', '');
        addOptional(p, 'isotherms', []);
        addOptional(p, 'plotIso', true);
        addOptional(p, 'fit', []);
        addOptional(p, 'canGetUfromgRind', []);
        addOptional(p, 'fileListOrgbyT', []);
        addOptional(p, 'hPvsRho', []);
        addOptional(p, 'hPvsV', []);
        addOptional(p, 'UVsT', true);
        addOptional(p, 'plotCv', true);
        addOptional(p, 'talk', false);
        parse(p, varargin{:});
        Results = p.Results;
        N = Results.N;
        Visible = Results.Visible;
        saveFig = Results.saveFig;
        fitprop = Results.fitprop;
        residuals = Results.residuals;
        fileNameEnd = Results.fileNameEnd;
        isotherms = Results.isotherms;
        fit = Results.fit;
        canGetUfromgRind = Results.canGetUfromgRind;
        fileListOrg = Results.fileListOrgbyT;
        hPvsRho = Results.hPvsRho;
        hPvsV = Results.hPvsV;
        UVsT = Results.UVsT;
        plotIso = Results.plotIso;
        talk = Results.talk;
        plotCv = Results.plotCv;
        
        if isempty(N) && isempty(isotherms) && isempty(fileListOrg)
            error(['provide number of particles,'...
                'for example: plotIso(''N'',625), or isotherm object, or fileListOrgbyT']);
        else
            fileList = dir(['N' num2str(N) 'T*mat']);
            fileList = {fileList.name};
        end
            
    
    if isempty(isotherms)
        if isempty(fileListOrg)
            fileListOrg = fileListOrgbyT(fileList);        
        end        
        
        [Niso, ~] = size(fileListOrg);
        
        for i = 1:Niso
            if i == 1
                isotherms = isotherm('datafileList',fileListOrg(i,:)','cutEquilirization',true);
            else
                isotherms(i) = isotherm('datafileList',fileListOrg(i,:)','cutEquilirization',true);
            end
        end
    else
        Niso = [];
    end
                    
                
    if isempty(Niso)
        Niso = length(isotherms);
    end
    Nrho = length(isotherms(1).rho);
    pressure = zeros(Niso,Nrho);
    rho = zeros(Niso,Nrho);
    
    for i = 1:Niso
        rho(i,:) = isotherms(i).rho;
        pressure(i,:) = isotherms(i).pressure;
        leg{1,i} = ['T = ' num2str(isotherms(i).T)];
        
    end
    
    if plotIso
        if isempty(hPvsRho)
            hPvsRho = figure('Visible',Visible);
        end
        hPvsRho = colorPlot(rho,pressure,'addLegend',leg,'figHandle',hPvsRho);
        title(['isotherms for N = ' num2str(N)]);
        xlabel('density, reduced units');
        ylabel('pressure, reduced units');

        if ~isempty(fitprop) && isempty(fit)
            for i = 1:length(fitprop)
                fit{1,i} = fitIso(isotherms,fitprop{i},'figureHandle',hPvsRho);
            end
        end

        if ~isempty(residuals)
            for i = 1:length(fit)
                if residuals{1,i}
                    for iso = 1:length(isotherms)
                        f = fit{1,i}(1,iso);
                        res = my_residuals(pressure(iso,:),f.P);
                        for j = 1:length(res)
                            text(rho(iso,j),f.P(j),num2str(res(j)));
                        end
                    end
                end
            end
        end

        if isempty(canGetUfromgRind)
            canGetUfromgRind = canGetUfromgR(isotherms,fit,0.1);
        end

        for i = 1:length(isotherms)
            plot(isotherms(i).rho(find(canGetUfromgRind(i,:))),...
                isotherms(i).pressure(find(canGetUfromgRind(i,:))),'mo');
        end

        if saveFig
            saveas(hPvsRho,['isotherms_N' num2str(N) 'P_rho' fileNameEnd '.fig']);
            saveas(hPvsRho,['isotherms_N' num2str(N) 'P_rho' fileNameEnd '.jpg']);
            if talk
                disp(['saved: isotherms_N' num2str(N) 'P_rho' fileNameEnd]);
            end
        end


        if isempty(hPvsRho)
            hPvsV = figure('Visible',Visible);
        end
        hPvsV = colorPlot(1./rho,pressure,'addLegend',leg,'figHandle',hPvsV);
        title(['isotherms for N = ' num2str(N)]);
        xlabel('volume, reduced units');
        ylabel('pressure, reduced units');
    
   
        if saveFig
            saveas(hPvsV,['isotherms_N' num2str(N) 'P_V' fileNameEnd '.fig']);
            saveas(hPvsV,['isotherms_N' num2str(N) 'P_V' fileNameEnd '.jpg']);
            
            if talk
                disp(['saved: isotherms_N' num2str(N) 'P_V' fileNameEnd]);
            end
        end
    end
    
    if UVsT
        U = zeros(Nrho,Niso);
        T = zeros(Nrho,Niso);
        N = num2str(isotherms(1).MC2DLJ(1).simulationParam.N);
        for i = 1:Niso
            for j = 1:Nrho
                U(j,i) = isotherms(i).MC2DLJ(j).data.meanUlrcEq;
                T(j,i) = isotherms(i).T;
                if i == 1
                    legUVsT{1,j} = ['\rho = ' my_num2str(isotherms(i).rho(j))];
                end
            end
        end
        
        colorPlot(T,U,'addLegend',legUVsT);
        title(['mean U Vs. T for N = ' N]);
        xlabel('Temperature, reduced units');
        ylabel('Energy, reduced units');

        if saveFig
            figName = ['meanUvsT_N' num2str(N) fileNameEnd];
            saveas(gcf, [figName '.fig']);
            saveas(gcf,[figName '.jpg']);
            
            if talk
                disp(['saved: ' figName]);
            end
        end


    end
    
    if plotCv
        cv = zeros(Nrho,Niso);
        T = zeros(Nrho,Niso);
        N = num2str(isotherms(1).MC2DLJ(1).simulationParam.N);
        for i = 1:Niso
            for j = 1:Nrho
                cv(j,i) = isotherms(i).MC2DLJ(j).data.cv;
                T(j,i) = isotherms(i).T;
                if i == 1
                    legUVsT{1,j} = ['\rho = ' my_num2str(isotherms(i).rho(j))];
                end
            end
        end
        
        colorPlot(T,cv,'addLegend',legUVsT);
        title(['cv Vs. T for N = ' N]);
        xlabel('Temperature, reduced units');
        ylabel('cv, reduced units');

        if saveFig
            figName = ['cvVsT_N' num2str(N) fileNameEnd];
            saveas(gcf, [figName '.fig']);
            saveas(gcf,[figName '.jpg']);
            
            if talk
                disp(['saved: ' figName]);
            end
        end
        
        cv = zeros(Niso, Nrho);
        rho = zeros(Niso, Nrho);
        N = num2str(isotherms(1).MC2DLJ(1).simulationParam.N);
        for i = 1:Niso
            legcvVsrho{1,j} = ['T = ' my_num2str(isotherms(i).T)];

            for j = 1:Nrho
                cv(i,j) = isotherms(i).MC2DLJ(j).data.cv;
                rho(i,j) = isotherms(i).rho;
            end
        end
        
        colorPlot(T,cv,'addLegend',legcvVsrho);
        title(['cv Vs. rho for N = ' N]);
        xlabel('rho, reduced units');
        ylabel('cv, reduced units');

        if saveFig
            figName = ['cvVsrho_N' num2str(N) fileNameEnd];
            saveas(gcf, [figName '.fig']);
            saveas(gcf,[figName '.jpg']);
            
            if talk
                disp(['saved: ' figName]);
            end
        end
    end
end

% functions used in this code
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

function ind = canGetUfromgR(isotherms,fit,maxres)
    
    ind = zeros(length(isotherms),length(isotherms(1).rho));

    for i = 1:length(isotherms)
        Plin = fit{1,1}(1,i).P;
        Pvir = fit{1,2}(1,i).P;
        Pdata = isotherms(i).pressure;
        
        reslin = abs((Plin-Pdata)./Pdata);
        resvir = abs((Pvir-Pdata)./Pdata);
        
        ind(i,:) = and((reslin >= maxres),(resvir <= maxres));
        
    end
end
