%save performence check
N = 600;
steps = 10000;
var = rand(2,N,steps);

% version 7.3 file (updating varible size)
[timeVarV7_3varSizeUpdate, stepsVarV7_3varSizeUpdate] =...
    saveVariable(var,'matFilev7_3varSizeUpdate.mat','version','-v7.3',...
    'updateVarSize',true,'saveBlocSize',1000);

% matFilev7_3 = matfile('matFilev7_3varSizeUpdate.mat');
% matFilev7_3 = matfile('matFilev7_3varSizeUpdate.mat','Writable',true);
% s = zeros(2,N,2);
% matFilev7_3.var = s;
% 
% timeVarV7_3varSizeUpdate = zeros(1,steps);
% for i = 1:steps
%     tic;
%     matFilev7_3.var(:,:,i) = var(:,:,i);
%     timeVarV7_3varSizeUpdate(1,i) = toc;
%     disp(i);
% end

% version 7.3 file (variable size set upfront)

[timeVarV7_3varNoSizeUpdate, stepsVarV7_3varNoSizeUpdate] =...
    saveVariable(var,'matFilev7_3varNoSizeUpdate.mat','version','-v7.3',...
    'updateVarSize',false,'saveBlocSize',1000);

% matFilev7_3 = matfile('matFilev7_3varNoSizeUpdate.mat');
% matFilev7_3 = matfile('matFilev7_3varNoSizeUpdate.mat','Writable',true);
% s = zeros(2,N,steps);
% matFilev7_3.var = s;
% 
% timeVarV7_3varNoSizeUpdate = zeros(1,steps);
% for i = 1:steps
%     tic;
%     matFilev7_3.var(:,:,i) = var(:,:,i);
%     timeVarV7_3varNoSizeUpdate(1,i) = toc;
%     disp(i);
% end



% version 7.3 file (savefast updating varible size)

% version 7.3 file (savefast size set upfront)

% version 7 file (updating varible size)

% version 7 file (variable size set upfront)

% version 7 file (savefast updating varible size)

% version 7 file (savefast size set upfront)

% version 6 file (updating varible size)

% version 6 file (variable size set upfront)

% version 6 file (savefast updating varible size)

% version 6 file (savefast size set upfront)

% fwrite

