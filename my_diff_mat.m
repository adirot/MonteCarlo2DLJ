function div = my_diff_mat(matx,maty)
% derivatives for every row

diffy = diff(maty,1,2); % diff between elements in every row
diffx = diff(matx,1,2);
div = diffy./diffx;

end