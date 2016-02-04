function [divx,divy] = my_diff_mat(matx,maty)
% derivatives for every row

diffy = diff(maty,1,2); % diff between elements in every row
diffx = diff(matx,1,2);
divy = diffy./diffx;

divx = (matx(:,2:end) + matx(:,1:(end-1)))/2;

end