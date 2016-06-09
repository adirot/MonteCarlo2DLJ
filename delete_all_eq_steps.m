%% delete first 4000 steps %%

start = 4000;

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
            
            if and(~(and(tind==1,and(rind==1,mind==1))),~(and(tind==1,and(rind==1,mind==2))))
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
                M = M.deleteFirstRuns(start);
                toc
            end
        end 
    end
end


