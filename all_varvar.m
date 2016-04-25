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
		allP = M.data.allPlrc;
		allU = M.data.allUlrc;
		RDFp = M.data.RDFhisto(:,31)';
		for k = 2:length(allP) stdP(k-1) = std(allP(1:k))/mean(allP(1:k)); end;
		for k = 2:length(allU) stdU(k-1) = std(allU(1:k))/mean(allU(1:k)); end;
		% find the first non zero element:
		ind = min(find(RDFp));
		for k = ind:length(allP) stdRDF(k-ind+1) = std(RDFp(1:k))/mean(RDFp(1:k)); end;
		
		for k = 2:length(stdP) stdstdP(k-1) = std(stdP(1:k))/mean(stdP(1:k)); end;
		for k = 2:length(stdP) stdstdU(k-1) = std(stdU(1:k))/mean(stdU(1:k)); end;
		for k = 2:length(stdRDF) stdstdRDF(k-1) = std(stdRDF(1:k))/mean(stdRDF(1:k)); end;
		
		allstdstdP{i,j} = stdstdP;
		allstdstdU{i,j} = stdstdU;
		allstdstdRDF{i,j} = stdstdRDF;

		stdstdendP(i,j) = stdstdP(end);
		stdstdendU(i,j) = stdstdU(end);
		stdstdendRDF(i,j) = stdstdRDF(end);
		
		figure; hold on;
		plot(stdstdP);
		plot(stdstdU,'r');
		plot(stdstdRDF,'k');
		
		title(['std of std U,P,RDF for T = ' num2str(T(i)) ' rho = ' num2str(rho(j))]);
		xlabel('recorded steps');
		ylabel('std of std');
		legend({'pressure', 'energy', 'RDF peak'});
		set(findall(gcf,'-property','FontSize'),'FontSize',22);
		
		saveas(gcf,['stdstdT' my_num2str(T(i)) 'rho' my_num2str(rho(j)) '.fig']);
		saveas(gcf,['stdstdT' my_num2str(T(i)) 'rho' my_num2str(rho(j)) '.jpg']);
		close all;
		disp(i);
		disp(j);
	end
end

save('stdstd.mat','stdstdP','stdstdU','stdstdRDF','stdstdendP','stdstdendU','stdstdendRDF');

		
