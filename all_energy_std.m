%% get all stds in energy %%


steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};

tind = 0;
rind = 0;
mind = 0;

for t = [0.45,0.6,0.8,1,1.5,2]
    tind = tind + 1;
    rind = 0;
    for r = [0.005,0.01,0.05,0.1]
        rind = rind + 1;
        mind = 0;
        for m = [3,4,5,6]
            mind = mind + 1; 
            
            disp(t); 
            disp(r);
            disp(m);
            disp('-----');
            
            list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
                '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
            end

            M = MC2DLJoutput(list{1,1});
tic;
 size(M.data.allUlrc)
            U10var{tind,rind,mind} = var(M.data.allUlrc(1,1:10));
            U50var{tind,rind,mind} = var(M.data.allUlrc(1,1:50));
            U100var{tind,rind,mind} = var(M.data.allUlrc(1,1:100));
            U1kvar{tind,rind,mind} = var(M.data.allUlrc(1,1:1000));
            U2kvar{tind,rind,mind} = var(M.data.allUlrc(1,1:2000));
            U6kvar{tind,rind,mind} = var(M.data.allUlrc(1,1:6000));
            U8kvar{tind,rind,mind} = var(M.data.allUlrc(1,1:8000));
        
            
            
            
       
%             plot(x,U12k{tind,rind,mind});            
%             hold on;
%             plot(x,U12k{tind,rind,mind}+...
%                 (U12k{tind,rind,mind}.*U12kstd{tind,rind,mind})/sqrt(12000),'-r');
%             
%             plot(x,U12k{tind,rind,mind}-...
%                 (U12k{tind,rind,mind}.*U12kstd{tind,rind,mind})/sqrt(12000),'-r');
%             
%             plot(x,U12kvar{tind,rind,mind},'r');
%             plot(x,U12kmax{tind,rind,mind},'r');
%             
%             saveas(gcf,['Ustdvarmax12kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax12kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
% 
%             
%             plot(x,U10{tind,rind,mind});            
%             hold on;
%             plot(x,U10{tind,rind,mind}+...
%                 (U10{tind,rind,mind}.*U10std{tind,rind,mind})/sqrt(10),'-r');
%             
%             plot(x,U10{tind,rind,mind}-...
%                 (U10{tind,rind,mind}.*U10std{tind,rind,mind})/sqrt(10),'-r');
%             
%             plot(x,U10var{tind,rind,mind},'r');
%             plot(x,U10max{tind,rind,mind},'r');
%             
%             saveas(gcf,['Ustdvarmax10T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax10T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,U50{tind,rind,mind});            
%             hold on;
%             plot(x,U50{tind,rind,mind}+...
%                 (U50{tind,rind,mind}.*U50std{tind,rind,mind})/sqrt(50),'-r');
%             
%             plot(x,U50{tind,rind,mind}-...
%                 (U50{tind,rind,mind}.*U50std{tind,rind,mind})/sqrt(50),'-r');
%             
%             
%             plot(x,U50var{tind,rind,mind},'r');
%             plot(x,U50max{tind,rind,mind},'r');
%             
%             saveas(gcf,['Ustdvarmax50T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax50T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,U100{tind,rind,mind});            
%             hold on;
%             plot(x,U100{tind,rind,mind}+...
%                 (U100{tind,rind,mind}.*U100std{tind,rind,mind})/sqrt(100),'-r');
%             
%             plot(x,U100{tind,rind,mind}-...
%                 (U100{tind,rind,mind}.*U100std{tind,rind,mind})/sqrt(100),'-r');
%             
%             
%             plot(x,U100var{tind,rind,mind},'r');
%             plot(x,U100max{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['Ustdvarmax100T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax100T' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,U1k{tind,rind,mind});            
%             hold on;
%             plot(x,U1k{tind,rind,mind}+...
%                 (U1k{tind,rind,mind}.*U1kstd{tind,rind,mind})/sqrt(1000),'-r');
%             
%             plot(x,U1k{tind,rind,mind}-...
%                 (U1k{tind,rind,mind}.*U1kstd{tind,rind,mind})/sqrt(1000),'-r');
%             
%             
%             plot(x,U1kvar{tind,rind,mind},'r');
%             plot(x,U1kmax{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['Ustdvarmax1kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax1kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,U2k{tind,rind,mind});            
%             hold on;
%             plot(x,U2k{tind,rind,mind}+...
%                 (U2k{tind,rind,mind}.*U2kstd{tind,rind,mind})/sqrt(2000),'-r');
%             
%             plot(x,U2k{tind,rind,mind}-...
%                 (U2k{tind,rind,mind}.*U2kstd{tind,rind,mind})/sqrt(2000),'-r');
%             
%             
%             plot(x,U2kvar{tind,rind,mind},'r');
%             plot(x,U2kmax{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['Ustdvarmax2kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax2kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             plot(x,U6k{tind,rind,mind});            
%             hold on;
%             plot(x,U6k{tind,rind,mind}+...
%                 (U6k{tind,rind,mind}.*U6kstd{tind,rind,mind})/sqrt(6000),'-r');
%             
%             plot(x,U6k{tind,rind,mind}-...
%                 (U6k{tind,rind,mind}.*U6kstd{tind,rind,mind})/sqrt(6000),'-r');
%             
%             
%             plot(x,U6kvar{tind,rind,mind},'r');
%             plot(x,U6kmax{tind,rind,mind},'r');
%             
%             
%             saveas(gcf,['Ustdvarmax6kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.fig']);
%             saveas(gcf,['Ustdvarmax6kT' my_num2str(t)...
%                 'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
%             close all;
%             
%             
%             close all;
        toc
        end 
    end
end

save('all_U_var.mat');
