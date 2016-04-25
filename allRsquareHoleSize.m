% get all rsquares and hole size

steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
start = 4000;

tind = 0;
rind = 0;
mind = 0;

rsquare10 = zeros(6,4,4);
rsquare50 = zeros(6,4,4);
rsquare100 = zeros(6,4,4);
rsquare1k = zeros(6,4,4);
rsquare2k = zeros(6,4,4);
rsquare6k = zeros(6,4,4);
rsquare12k = zeros(6,4,4);

umin10 = zeros(6,4,4);
umin50 = zeros(6,4,4);
umin100 = zeros(6,4,4);
umin1k = zeros(6,4,4);
umin2k = zeros(6,4,4);
umin6k = zeros(6,4,4);
umin12k = zeros(6,4,4);

for t = [0.45,0.6,0.8,1,1.5,2]
    tind = tind + 1;
	rind = 0;
	mind = 0;
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

            if and(and(t == 0.45,r == 0.005), m==3)
                x = M.data.RDFbins;
            end

            RDF10 = mean(M.data.RDFhisto(start:(start+10),1:300));
            RDF50 = mean(M.data.RDFhisto(start:(start+50),1:300));
            RDF100 = mean(M.data.RDFhisto(start:(start+100),1:300));
            RDF1k = mean(M.data.RDFhisto(start:(start+1000),1:300));
            RDF2k = mean(M.data.RDFhisto(start:(start+2000),1:300));
            RDF6k = mean(M.data.RDFhisto(start:(start+6000),1:300));
            RDF12k = mean(M.data.RDFhisto(start:(start+12000),1:300));

            [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF10), t, r,m,steps,'plotFig', false);
            rsquare10(tind,rind,mind) = gof.rsquare;
            umin10(tind,rind,mind) =...
                4*((mfit/12)^(-12/(mfit-12))-(mfit/12)^(-mfit/(mfit-12)));
            

            [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF50), t, r,m,steps,'plotFig', false);
            rsquare50(tind,rind,mind) = gof.rsquare;
            umin50(tind,rind,mind) =...
                4*((mfit/12)^(-12/(mfit-12))-(mfit/12)^(-mfit/(mfit-12)));
               
            [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF100), t, r,m,steps,'plotFig', false);
            rsquare100(tind,rind,mind) = gof.rsquare;
            umin100(tind,rind,mind) =...
                4*((mfit/12)^(-12/(mfit-12))-(mfit/12)^(-mfit/(mfit-12)));
            
            [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF1k), t, r,m,steps,'plotFig', false);
            rsquare1k(tind,rind,mind) = gof.rsquare;
            umin1k(tind,rind,mind) =...
                4*((mfit/12)^(-12/(mfit-12))-(mfit/12)^(-mfit/(mfit-12)));
            
            
            [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF2k), t, r,m,steps,'plotFig', false);
            rsquare2k(tind,rind,mind) = gof.rsquare;
            umin2k(tind,rind,mind) =...
                4*((mfit/12)^(-12/(mfit-12))-(mfit/12)^(-mfit/(mfit-12)));
                        [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF6k), t, r,m,steps,'plotFig', false);
            rsquare6k(tind,rind,mind) = gof.rsquare;

            [fitresult, mfit, mError, gof] =...
                createFitRDF(x, -log(RDF12k), t, r,m,steps,'plotFig', false);
            rsquare12k(tind,rind,mind) = gof.rsquare;   
            umin12k(tind,rind,mind) =...
                4*((mfit/12)^(-12/(mfit-12))-(mfit/12)^(-mfit/(mfit-12)));
                        
        end
    end
end

save('allRsquaredHolesize.mat','umin10','umin50','umin100','umin1k'...
    ,'umin2k','umin6k','umin12k'...
    ,'rsquare10','rsquare50','rsquare100','rsquare1k',...
    'rsquare2k','rsquare6k','rsquare12k');
