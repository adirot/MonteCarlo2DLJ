%% get all U from RDF %%

tind = 0;

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
            
            if ~and(and(tind==1,rind==1),mind==1)
                list = dir(['*T' my_num2str(t) 'rho' my_num2str(r) '*_m' num2str(m)...
                    '*mat']);
                list = {list.name}';

                if length(list) > 1
                    disp('more than one file!');
                end

                M = MC2DLJoutput(list{1,1});
                allmeanU{tind,rind,mind} = mean(M.data.allU);
                
                if existInMatfile(M.data, 'UfromRDFinStep');
                    allmeanUfromRDF{tind,rind,mind} = mean(M.data.UfromRDFinStep);
                    allUfromRDF{tind,rind,mind} = M.data.UfromRDFinStep;
                else
                    [M, ~] = M.getVarUfromRDF('use_m_n_T',...
                        [M.simulationParam.m, 12, M.simulationParam.T]);
                    allUfromRDF{tind,rind,mind} = M.data.UfromRDFinStep;
                end
                
                [~, lenStepInd] = size(M.data.stepInd);
                [~, lenAllU] = size(M.data.allU);
                [~, lenallUfromRDF] = size(allUfromRDF{tind,rind,mind});
                
                if lenStepInd > lenAllU
                    plot(M.data.stepInd(1,(lenStepInd-lenAllU+1):lenStepInd)/M.simulationParam.N,...
                        M.data.allU);
                else
                    if lenStepInd < lenAllU
                        plot(M.data.stepInd/M.simulationParam.N,...
                            M.data.allU(1,(lenAllU-lenStepInd+1):lenAllU));
                    else
                        plot(M.data.stepInd/M.simulationParam.N,...
                            M.data.allU);
                    end
                end
                
                hold on;
                
                if lenStepInd > lenallUfromRDF
                    plot(M.data.stepInd(1,(lenStepInd-lenallUfromRDF+1):lenStepInd)/M.simulationParam.N,...
                        allUfromRDF{tind,rind,mind},'r');
                else
                    if lenStepInd < lenallUfromRDF
                        plot(M.data.stepInd/M.simulationParam.N,...
                            allUfromRDF{tind,rind,mind}(1,(lenallUfromRDF-lenStepInd+1):lenallUfromRDF),'r');
                    else
                        plot(M.data.stepInd/M.simulationParam.N,...
                            allUfromRDF{tind,rind,mind},'r');
                    end
                end
                
                xlabel('Sweeps');
                ylabel('U');
                title('U from simulation compared with U from integrating over RDF');
                legend({'U from RDF', 'U from simulation'});
                
                saveas(gcf,['UfromRDF_T' my_num2str(t)...
                    'rho' my_num2str(r) 'm' num2str(m) '.fig']);
                saveas(gcf,['UfromRDF_T' my_num2str(t)...
                    'rho' my_num2str(r) 'm' num2str(m) '.jpg']);
                close all;
            end
        end 
    end
end

save('all_U_from_RDF.mat');
