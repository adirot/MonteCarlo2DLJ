clear all;
T = [0.45, 0.6, 0.8, 2,3,4,5,6,7,8];
rho = [0.01,0.1,0.3];

for i = 1:length(T)
		for j = 1:length(rho)

		list = dir(['N625T' my_num2str(T(i)) 'rho' my_num2str(rho(j)) 'i*off10*mat']);
		list = {list.name}';
		if length(list) > 1
			disp(['problem at T = ' num2str(T(i)) ' rho = ' num2str(rho(j))]);
		end
		M = MC2DLJoutput(list{1,1});
		if i==1
			if j == 1
				bins = M.data.RDFbins;
			end
		end
		M = MC2DLJoutput(list{1,1});
		
		first = 4000;
	
		RDF10{i,j} = mean(M.data.RDFhisto(first:(first+10),1:300));
		RDF100{i,j} = mean(M.data.RDFhisto(first:(first+100),1:300));
		RDF1k{i,j} = mean(M.data.RDFhisto(first:(first+1000),1:300));
		RDF2k{i,j} = mean(M.data.RDFhisto(first:(first+2000),1:300));
		RDF6k{i,j} = mean(M.data.RDFhisto(first:(first+6000),1:300));
		RDF12k{i,j} = mean(M.data.RDFhisto(first:(first+12000),1:300));

		disp(i);		
		disp(j);
	end
end

save('a.mat','RDF10','RDF100','RDF1k','RDF2k','RDF6k','RDF12k','bins');

		
