
indm = 0;

for m = [3,4,5,6]
	indm = indm + 1;
	
	for T = [0.45 0.6 0.8 1 1.5 2]
		if T == 0.45
			list = dir(['N*T' my_num2str(T) 'r*m' num2str(m) '*mat']);
			list = {list.name};
		end
		list2 = dir(['N*T' my_num2str(T) 'r*m' num2str(m) '*mat']);
		list2 = {list2.name};
		list = [list; list2];
	end
	
	for s = [10 50 100 1000 2000 6000 12000]
		[isotherms,fit,res,canGetUfromgRind,~,~,P,U,T,Z,Zx] = ...
			plotIso('fileListOrg',list,'N',625,'fitprop',{'plotLin', 'plotVirialExp'},...
			'residuals',{true, true},'talk',true,'plotCv',false,'endBy',4000+s);
		allres{1,indm,s} = res;
	end

end
save('allresiduals.mat','allres');
