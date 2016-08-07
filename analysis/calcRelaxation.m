function [steps,phi] = calcRelaxation(MC)
% the MC object should be so that every object is a different MC run

P = MC(1).data.allPlrc;
steps = MC(1).data
for i = 2:length(MC)
    P = P + MC(i).data.allPlrc;
end

Pmean = P/length(MC);
PfinMean = P(end)/length(MC);
PinitMean = P(1)/length(MC);

phi = (Pmean - PfinMean)/(PinitMean - PfinMean);
end
        
%i = 1; for t = [1 9 17 25 33] for j = t:(t+7) tau = trapz(steps(j,:),phi(j,:)); leg{1,j} = ['\rho = ' num2str(rho(j)) ' \tau = ' num2str(tau)];end; colorPlot(steps(t:(t+7),:),phi(t:(t+7),:),'addLegend',leg(1,t:(t+7))); name = ['relaxPN625T' my_num2str(T(i))];i = i+1; xlabel('steps'); ylabel('\phi - relaxation parameter'); title(name); saveas(gcf,[name '.fig']); saveas(gcf,[name '.jpg']);disp(t); close all; end;
