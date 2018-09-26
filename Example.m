clear; clc;

%% design network
Net1 = BrainNetSim(3,300,.96);

% define net dynamics
Net1 = Net1.AddNodeFreqs([1 2 3],{[8],[5 10],[4 15]});

Net1 = Net1.AddConnection([1 2],'Type','bandpass','LF',5,'HF',20,'Gain',1);

Net1 = Net1.GenerateARMatrix;
NS = 20000; % length of simulation: time points
[Net1,TS] = Net1.Realization(NS);


%%
figure,

for node = 1:Net1.NodeNum
    % 
    subplot(Net1.NodeNum,2,(node-1)*2+1),plot(TS(node,1:Net1.SF*2));
    set(gca,'xtick',0:Net1.SF:Net1.SF*2,'xticklabel',0:2)
    xlabel('time (s)')
    ylabel(['Node' num2str(node)],'fontweight','bold','fontsize',12);
    title('Temporal dynamic');
     
    subplot(Net1.NodeNum,2,(node-1)*2+2);
    [Z,f] = pwelch(TS(node,:),Net1.SF,[],[],Net1.SF);
    plot(f,Z,'linewidth',2);title('PSD');xlabel('Frequency(Hz)');xlim([0 50])
    %legend(['Signal #' num2str(Sig)])
end

