function res = my_residuals(data,fit)
    res = (data-fit)./data;
end