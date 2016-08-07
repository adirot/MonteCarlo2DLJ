function my_var = my_var_from_std(my_std)

for tind = 1:6
    for rind = 1:4
        for mind = 1:4
            my_var{tind,rind,mind} = my_std{tind,rind,mind}.^2;
        end
    end
end

end