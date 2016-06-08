% get RDFs with different bin sizes

start = 4000;

tind = 0;
rind = 0;
mind = 0;

timeel = [];

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
	    tic;
            
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

	    [M, bins50, RDF50] = M.calcRDF(10,50,'save2data',false,'startFrom',start,'talk',false);		
	    disp('50');
	    [M, bins100, RDF100] = M.calcRDF(10,100,'save2data',false,'startFrom',start,'talk',false);		
	    disp('100');
	    [M, bins150, RDF150] = M.calcRDF(10,150,'save2data',false,'startFrom',start,'talk',false);		
 	    disp('150');
	    [M, bins200, RDF200] = M.calcRDF(10,200,'save2data',false,'startFrom',start,'talk',false);		
	    disp('200');
	    [M, bins250, RDF250] = M.calcRDF(10,250,'save2data',false,'startFrom',start,'talk',false);		
	    disp('250');
	    [M, bins350, RDF350] = M.calcRDF(10,350,'save2data',false,'startFrom',start,'talk',false);		
	    disp('350');
	    [M, bins400, RDF400] = M.calcRDF(10,400,'save2data',false,'startFrom',start,'talk',false);		
	    disp('400');

	    bins300 = M.data.RDFbins;
	    RDF300 = M.data.RDFhisto;		
            RDF300 = RDF300(start:(start+12000),:);
	    disp('300');
	    
            x = [...
		padarray(bins50',400-50,'post')';...
		padarray(bins100',400-100,'post')';...
		padarray(bins150',400-150,'post')';...
		padarray(bins200',400-200,'post')';...
		padarray(bins250',400-250,'post')';...
		padarray(bins300',400-300,'post')';...
		padarray(bins350',400-350,'post')';...
		bins400];
            y = [...
		padarray(RDF50',400-50,'post')';...
		padarray(RDF100',400-100,'post')';...
		padarray(RDF150',400-150,'post')';...
		padarray(RDF200',400-200,'post')';...
		padarray(RDF250',400-250,'post')';...
		padarray(RDF300',400-300,'post')';...
		padarray(RDF350',400-350,'post')';...
		RDF400];
		
	    colorPlot(x,y,'addLegend',{'50','100','150','200','250','300','350','400'},'length2plot',[50,100,150,200,250,300,350,400]);
	    title(['RDF diff number of bins 50-400. T = ' num2str(t) ' \rho = ' num2str(r) ' m = ' num2str(m)]);
	    xlabel('r');
	    ylabel('g(r)');
	    saveas(gcf,['RDFdiffBins50_400T' my_num2str(t) 'rho' my_num2str(r) 'm' my_num2str(m) '.fig']);             
	    saveas(gcf,['RDFdiffBins50_400T' my_num2str(t) 'rho' my_num2str(r) 'm' my_num2str(m) '.jpg']);             
	    close all;
	    disp('ploted and saved');
	    toc
	    timeel = [timeel toc];
	    save('time_run_diff_bin_RDF.mat','timeel');
        end
    end
end

