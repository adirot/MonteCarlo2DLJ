
% fit all matfiles to the potantail (logRDF) for different steps

steps = {'10';'50'; '100'; '1k'; '2k'; '6k'; '12k'};
list = dir('N*mat');
list = {list.name}';
start = 4000;

tind = 0;
rind = 0;
mind = 0;
nind = 0;

for t = [0.45,0.6,0.8,1,1.5,2]
    tind = tind + 1;
    rind = 0;
    for r = [0.005,0.01,0.05,0.1]
        rind = rind + 1;
        mind = 0;
        for m = [3,4,5,6]
            mind = mind + 1; 
	nind  = 0;
        for n = 7:16
            nind = nind + 1;

            disp(t); 
            disp(r);
            disp(m);
            disp(n);
            disp('-----');
            
            list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
                '*mat']);
            list = {list.name}';

            if length(list) > 1
                disp('more than one file!');
            end

            M = MC2DLJoutput(list{1,1});

            if and(and(and(t == 0.45,r == 0.005), m==3),n==7)
                x = M.data.RDFbins;
            end
            xs = [x;x;x;x;x;x;x;];
            
            tic;     
        
            RDF10 = mean(M.data.RDFhisto(start:(start+10),1:300));
            RDF50 = mean(M.data.RDFhisto(start:(start+50),1:300));
            RDF100 = mean(M.data.RDFhisto(start:(start+100),1:300));
            RDF1k = mean(M.data.RDFhisto(start:(start+1000),1:300));
            RDF2k = mean(M.data.RDFhisto(start:(start+2000),1:300));
            RDF6k = mean(M.data.RDFhisto(start:(start+6000),1:300));
            RDF12k = mean(M.data.RDFhisto(start:(start+12000),1:300));
        
            ys = [-log(RDF10);-log(RDF50);-log(RDF100);-log(RDF1k);...
                -log(RDF2k);-log(RDF6k);-log(RDF12k)];

            % with n set to 12
            %[fitresult,mfit,merror, nfit, nerror, Tfit, Terror, gof] = createFitRDF(xs, ys, t, r,m,steps);
        
            % with free n
%             [fitresult,mfit,merror,nfit,nerror,Tfit,Terror, gof] =...
%                 createFitRDF(xs, ys, t, r,m,steps,'freeTandn',true);
 
            % with free T, bound n
            % [fitresult,mfit,merror,nfit,nerror,Tfit,Terror, gof] =...
            %     createFitRDF(xs, ys, t, r,m,steps,'freeTnbound',true);

            % with free T, n set
            [fitresult,mfit,merror,nfit,nerror,Tfit,Terror, gof] =...
                 createFitRDF(xs, ys, t, r,m,steps,'freeTnset',true,'nset',n);
            

            allfittedm10{tind,rind,mind,nind} = mfit(1);
            allfittedm50{tind,rind,mind,nind} = mfit(2);
            allfittedm100{tind,rind,mind,nind} = mfit(3);
            allfittedm1k{tind,rind,mind,nind} = mfit(4);
            allfittedm2k{tind,rind,mind,nind} = mfit(5);
            allfittedm6k{tind,rind,mind,nind} = mfit(6);
            allfittedm12k{tind,rind,mind,nind} = mfit(7);
            
            allmerror10{tind,rind,mind,nind} = merror(1);
            allmerror50{tind,rind,mind,nind} = merror(2);
            allmerror100{tind,rind,mind,nind} = merror(3);
            allmerror1k{tind,rind,mind,nind} = merror(4);
            allmerror2k{tind,rind,mind,nind} = merror(5);
            allmerror6k{tind,rind,mind,nind} = merror(6);
            allmerror12k{tind,rind,mind,nind} = merror(7);

            %allfittedn10{tind,rind,mind} = nfit(1);
            %allfittedn50{tind,rind,mind} = nfit(2);
            %allfittedn100{tind,rind,mind} = nfit(3);
            %allfittedn1k{tind,rind,mind} = nfit(4);
            %allfittedn2k{tind,rind,mind} = nfit(5);
            %allfittedn6k{tind,rind,mind} = nfit(6);
            %allfittedn12k{tind,rind,mind} = nfit(7);
            
            %allnerror10{tind,rind,mind} = nerror(1);
            %allnerror50{tind,rind,mind} = nerror(2);
            %allnerror100{tind,rind,mind} = nerror(3);
            %allnerror1k{tind,rind,mind} = nerror(4);
            %allnerror2k{tind,rind,mind} = nerror(5);
            %allnerror6k{tind,rind,mind} = nerror(6);
            %allnerror12k{tind,rind,mind} = nerror(7);

            allfittedT10{tind,rind,mind,nind} = Tfit(1);
            allfittedT50{tind,rind,mind,nind} = Tfit(2);
            allfittedT100{tind,rind,mind,nind} = Tfit(3);
            allfittedT1k{tind,rind,mind,nind} = Tfit(4);
            allfittedT2k{tind,rind,mind,nind} = Tfit(5);
            allfittedT6k{tind,rind,mind,nind} = Tfit(6);
            allfittedT12k{tind,rind,mind,nind} = Tfit(7);
            
            allTerror10{tind,rind,mind,nind} = Terror(1);
            allTerror50{tind,rind,mind,nind} = Terror(2);
            allTerror100{tind,rind,mind,nind} = Terror(3);
            allTerror1k{tind,rind,mind,nind} = Terror(4);
            allTerror2k{tind,rind,mind,nind} = Terror(5);
            allTerror6k{tind,rind,mind,nind} = Terror(6);
            allTerror12k{tind,rind,mind,nind} = Terror(7);

            saveas(gcf,['fitRDFT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) 'freeTnset' num2str(n) '_steps.fig']);
            saveas(gcf,['fitRDFT' my_num2str(t)...
                'rho' my_num2str(r) 'm' num2str(m) 'freenTnset' num2str(n) '_steps.jpg']);
            close all;
       
        toc
        end
        end 
    end
end

save('all_fittedm_freeTnset7_16.mat');
