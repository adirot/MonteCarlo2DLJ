    function MC2DLJs = getMC2DLJs(fileListOrgbyT)
        [Niso, Nrho] = size(fileListOrgbyT);
        
        for i = 1:Niso
            for j = 1:Nrho
                if i == 1 && j == 1
                    MC2DLJs = MC2DLJoutput(fileListOrgbyT{1,1});
                else
			if exist(fileListOrgbyT{i,j},'file') == 2
                    		MC2DLJs(i,j) = MC2DLJoutput(fileListOrgbyT{i,j});
			end
                end
            end
        end
    end