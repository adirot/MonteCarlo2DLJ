% std for m dirived for 10\50\100 steps.

steps = {'10';'50'; '100'; '1k';};
start = 4000;

mfit10 = [];
mfit50 = [];
mfit100 = [];
mfit1k = [];

allstds10 = zeros(6,4,4);
allstds50 = zeros(6,4,4);
allstds100 = zeros(6,4,4);
allstds1k = zeros(6,4,4);

allmeans10 = zeros(6,4,4);
allmeans50 = zeros(6,4,4);
allmeans100 = zeros(6,4,4);
allmeans1k = zeros(6,4,4);


tind = 0;
rind = 0;
mind = 0;

for t = [0.45,0.6,0.8,1,1.5,2]
    tind = tind + 1;
    for r = [0.005,0.01,0.05,0.1]
        rind = rind + 1;
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

            if and(and(t == 0.45,r == 0.005), m==3)
                x = M.data.RDFbins;
            end

            for j = 1:10
            disp(j);
            tic
                rndind = unidrnd(12000,1,10)+start;
                RDF10 = [];
                for k = 1:length(rndind)
                    RDF10 = [RDF10; M.data.RDFhisto(rndind(k),1:300)];
                end
                RDF10 = mean(RDF10);
                
                rndind = unidrnd(12000,1,50)+start;
                RDF50 = [];
                for k = 1:length(rndind)
                    RDF50 = [RDF50; M.data.RDFhisto(rndind(k),1:300)];
                end
                RDF50 = mean(RDF50);
                
                rndind = unidrnd(12000,1,100)+start;
                RDF100 = [];
                for k = 1:length(rndind)
                    RDF100 = [RDF100; M.data.RDFhisto(rndind(k),1:300)];
                end
                RDF100 = mean(RDF100);
                
                rndind = unidrnd(12000,1,1000)+start;
                RDF1k = [];
                for k = 1:length(rndind)
                    RDF1k = [RDF1k; M.data.RDFhisto(rndind(k),1:300)];
                end
                RDF1k = mean(RDF1k);
            
                [fitresult, mfit, mError, gof] =...
                    createFitRDF(x, -log(RDF10), t, r,m,steps,'plotFig', false);
                mfit10 = [mfit10 mfit];

                [fitresult, mfit, mError, gof] =...
                    createFitRDF(x, -log(RDF50), t, r,m,steps,'plotFig', false);
                mfit50 = [mfit50 mfit];
               
                [fitresult, mfit, mError, gof] =...
                    createFitRDF(x, -log(RDF100), t, r,m,steps,'plotFig', false);
                mfit100 = [mfit100 mfit];
                
                [fitresult, mfit, mError, gof] =...
                    createFitRDF(x, -log(RDF1k), t, r,m,steps,'plotFig', false);
                mfit1k = [mfit1k mfit];
                toc             
            end
            
            allmean10(tind,rind,mind) = mean(mfit10);
            allstds10(tind,rind,mind) = std(mfit10);
            allmean50(tind,rind,mind) = mean(mfit50);
            allstds50(tind,rind,mind) = std(mfit50);
            allmean100(tind,rind,mind) = mean(mfit100);
            allstds100(tind,rind,mind) = std(mfit100);
            allmean1k(tind,rind,mind) = mean(mfit1k);
            allstds1k(tind,rind,mind) = std(mfit1k);
               
            
            % plot

            [x,y] = hist(mfit10);
            plot(x,y);
            hold on;
            [x,y] = hist(mfit50);
            plot(x,y,'r');
            [x,y] = hist(mfit100);
            plot(x,y,'k');
            [x,y] = hist(mfit1k);
            plot(x,y,'m');
            legend(steps);
            title(['m from fit to -log(RDF), for different number of steps. T = '...
                my_num2str(t) ' \rho = ' num2str(r) ' m = ' num2str(m)]);
            xlabel('m');
            ylabel('count');

            saveas(gcf,['m_from_fitRDF' my_num2str(t)...
                'r' my_num2str(r) 'm' num2str(m) 'diff_steps.fig']);
            saveas(gcf,['m_from_fitRDF' my_num2str(t)...
                'r' my_num2str(r) 'm' num2str(m) 'diff_steps.jpg']);
            close all;
            
            
            mfit10 = [];
            mfit50 = [];
            mfit100 = [];
            mfit1k = [];
        end
    end
end

save('stdsform_diff_steps.mat','allstds10','allstds50','allstds100','allstds1k'...
    ,'allmeans10','allmeans50','allmeans100','allmeans1k');
