% T_div_rho = [];
% RDFvarforplot12k = [];

%     for mind = 1:4
%         i = 0;
%         for tind = [0.45 0.6 0.8 1 1.5 2]
%             for rind =[0.005 0.01 0.05 0.1]
%                 i = i +1;
%                 Tsquared_div_rho{1,mind,i} = tind^2/rind;
%                 T_div_rho{1,mind,i} = tind/rind;
%             end
%         end
%     end

%     for mind = 1:4
%         i = 0;
%         for tind = [0.45 0.6 0.8 1 1.5 2]
%             for rind =[0.005 0.01 0.05 0.1]
%                 i = i +1;
%                 T_rho{1,mind,i} = tind*rind;
%                 
%             end
%         end
%     end

% 
% for j = 1:length(xind)
%     for mind = 1:4
%         i = 0;
%         for tind = 1:6
%             for rind = 1:4
%                 i = i +1;
%                 RDFvarforplot12k{j,mind,i} = RDF12kstd{tind,rind,mind}(1,xind(j));
%             end
%         end
%     end
% end
% 
% for j = 1:length(xind)
%     figure;
%     hold on; 
%     plot([T_div_rho{1,1,:}],[RDFvarforplot12k{j,1,:}],'.')
%     plot([T_div_rho{1,2,:}],[RDFvarforplot12k{j,2,:}],'.r')
%     plot([T_div_rho{1,3,:}],[RDFvarforplot12k{j,3,:}],'.k')
%     plot([T_div_rho{1,4,:}],[RDFvarforplot12k{j,4,:}],'.g')
%     legend('m = 3','m = 4','m = 5','m = 6');
%     xlabel('T/rho');
%     ylabel('RDF var');
%     title(['RDF var vs. T/rho different m, x = ' num2str(r(xind(j)))]);
% end
% 
% 
% for j = 1:length(xind)
%     figure;
%     hold on; 
%     plot([Tsquared_div_rho{1,1,:}],[RDFvarforplot12k{j,1,:}],'.')
%     plot([Tsquared_div_rho{1,2,:}],[RDFvarforplot12k{j,2,:}],'.r')
%     plot([Tsquared_div_rho{1,3,:}],[RDFvarforplot12k{j,3,:}],'.k')
%     plot([Tsquared_div_rho{1,4,:}],[RDFvarforplot12k{j,4,:}],'.g')
%     legend('m = 3','m = 4','m = 5','m = 6');
%     xlabel('T^2/rho');
%     ylabel('RDF var');
%     title(['RDF var vs. T^2/rho different m, x = ' num2str(r(xind(j)))]);
% end

% 
% for j = 1:length(xind)
%     figure;
%     hold on; 
%     plot([T_rho{1,1,:}],[RDFvarforplot12k{j,1,:}],'.')
%     plot([T_rho{1,2,:}],[RDFvarforplot12k{j,2,:}],'.r')
%     plot([T_rho{1,3,:}],[RDFvarforplot12k{j,3,:}],'.k')
%     plot([T_rho{1,4,:}],[RDFvarforplot12k{j,4,:}],'.g')
%     legend('m = 3','m = 4','m = 5','m = 6');
%     xlabel('T*rho');
%     ylabel('RDF var');
%     title(['RDF var vs. T*rho different m, x = ' num2str(r(xind(j)))]);
% end

%   for mind = 1:4
%         i = 0;
%         for tind = [0.45 0.6 0.8 1 1.5 2]
%             for rind =[0.005 0.01 0.05 0.1]
%                 i = i +1;
%                 T{1,mind,i} = tind;
%                 
%             end
%         end
%   end
%     
% for j = 1:length(xind)
%     figure;
%     hold on; 
%     plot([T{1,1,:}],[RDFvarforplot12k{j,1,:}],'.')
%     plot([T{1,2,:}],[RDFvarforplot12k{j,2,:}],'.r')
%     plot([T{1,3,:}],[RDFvarforplot12k{j,3,:}],'.k')
%     plot([T{1,4,:}],[RDFvarforplot12k{j,4,:}],'.g')
%     legend('m = 3','m = 4','m = 5','m = 6');
%     xlabel('T');
%     ylabel('RDF var');
%     title(['RDF var vs. T different m, x = ' num2str(r(xind(j)))]);
% end
% 
%  for mind = 1:4
%         i = 0;
%         for tind = [0.45 0.6 0.8 1 1.5 2]
%             for rind =[0.005 0.01 0.05 0.1]
%                 i = i +1;
%                 T{1,mind,i} = tind;
%                 
%             end
%         end
%  end
%     
%   


col = {'b', 'r','k','g'};
figure;
T = [0.45 0.6 0.8 1 1.5 2];
rho = [0.005 0.01 0.05 0.1];
for rind = 1:4
    for mind = 1:4
        for tind = 1:6
            
            hold on; 
            plot((T.^2),[U8kvar{:,rind,mind}],'.','color',col{1,rind});
            xlabel('T^2');
            ylabel('U var');
            title('U var vs. T^2 different m');
        end
    end 
end