clear; clc;

%% design network
Net1 = BrainNetSim(3);

% define net dynamics
Net1 = Net1.AddNodeFreqs([1 2 3],{[10],[5 20],[4]});

Net1 = Net1.AddConnection([1 2],'Type','bandpass','LF',5,'HF',20,'Gain',1);

Net1 = Net1.GenerateARMatrix;
NS = 10000; % length of simulation: time points
[Net1,TS] = Net1.Realization(NS);


%%
figure,

for node = 1:Net1.NodeNum
    % 
    subplot(Net1.NodeNum,2,(node-1)*2+1),plot(TS(node,1:Net1.SF*10))
    title(['Dynamic of Node' num2str(node)],'fontweight','bold','fontsize',12)
    subplot(Net1.NodeNum,2,(node-1)*2+2);
    [Z,f] = pwelch(TS(node,:),100,[],[],Net1.SF);
    plot(f,Z,'linewidth',2);ylabel('PSD');xlabel('Frequency(Hz)');xlim([0 Net1.SF/4])
    %legend(['Signal #' num2str(Sig)])
end

