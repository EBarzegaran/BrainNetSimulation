classdef BrainNetSim
    
    properties
        NodeNum
        Nodes % a strcut array of nodes (NodeNum x 1), storing nodes information including
              % node.Type, Node.Freq, Node.TS %timeSeries
        Connections % a strut matrix of Connections information (NodeNum x NodeNum)including
              % Edge.FiltType, Edge.FiltOrder, Edge.FiltLF,
              % Edge.FiltHF,Edge.Gain
        SF    % sampling frequency for realization of time series
        
        ARMatrix % autoregressive matrix
    end
    
    properties (Dependent)
       
    end
    
    methods
        
        %% network initialization
        function obj = BrainNetSim(NodeNum,SF,rlevel)
            % INPUT: 
                    % NodeNum: number of the nodes in the network
                    % SF: sampling frequency of the system
                    % rLevel: the amplitude of the node poles, value
                    % between 0 and 1, the closer to 1, the more
                    % oscillatory, default value .95
            obj.NodeNum = NodeNum;
            if ~exist('rlevel','var') || isempty(rlevel)
                rlevel = 0.95; 
            end
            obj.Nodes = struct('Type',cell(NodeNum,1),'Freq',cell(NodeNum,1),'TS',cell(NodeNum,1),'r',num2cell(ones(NodeNum,1)*rlevel));
            obj.Connections = struct('FiltType',cell(NodeNum),'FiltOrder',cell(NodeNum),'LF',cell(NodeNum),'HF',cell(NodeNum),'Gain',cell(NodeNum));
            if ~exist('SF','var') || isempty(SF)
                SF = 100;
            end
            obj.SF = SF;
        end
        
        %% assign internal frequencies
        function obj = AddNodeFreqs (obj,nodes, Freqs)
            % syntax NetStrct.AddNodeFreqs(Freqs)
                % Freqs: is a cell array of internal frequencies of the
                        % nodes, each of the cell should contain one or two
                        % frequencies
            %[obj.Nodes.Freq] = deal(Freqs{:});%, 4}; 
            for i = 1:numel(nodes)
                if ~isempty(Freqs{i})
                    obj.Nodes(nodes(i)).Type = 'AR';
                    obj.Nodes(nodes(i)).Freq = Freqs{i};
                end
            end
        end
        
        %% add Connections to network
        function obj = AddConnection(obj, nodes,varargin)
            % array of 2x1, indicating the node connections: ATTENTION
            % order is important, it indicates the direction of connection
            %
            opt = ParseArgs(varargin,...
                'Type'      ,'',...
                'LF'        ,[],...
                'HF'        ,[],...
                'Gain'      ,1, ...
                'Order'     ,6 ...
                );
            
            % check if the connection exist in the network
            if numel(nodes)~=2 || sum(~ismember(nodes, 1:obj.NodeNum))
                error ('nodes should be a 2x1 vector indicating the node number of the network');
            end
            %
            if ~ismember(opt.Type,{'low','high','bandpass'})
                error ('No such filter is defined');
            else
                obj.Connections(nodes(1),nodes(2)).FiltType = opt.Type;
            end
            obj.Connections(nodes(1),nodes(2)).LF = opt.LF;
            obj.Connections(nodes(1),nodes(2)).HF = opt.HF;
            obj.Connections(nodes(1),nodes(2)).Gain = opt.Gain; 
            obj.Connections(nodes(1),nodes(2)).FiltOrder = opt.Order; 
        end
        
        %% AR calculation for linear models only
        function obj = GenerateARMatrix(obj)
            
            % First find the parameters including r(s) and noise level         
            order = {obj.Connections.FiltOrder};order = cat(3,order{:});
            order = max(max(order)+1,4);

            if ~isempty(order) 
            else
                order = 4;
            end
            A = zeros(obj.NodeNum,obj.NodeNum,order);
            
            for n = 1: obj.NodeNum
                f = obj.Nodes(n).Freq/obj.SF; % first frequency
                [r,f]=parameter_search([2*pi*f],obj.Nodes(n).r,1);
                [x,y] = pol2cart(f,r);Poles = x+1i*y;
                [~,a] = zp2tf(0,[Poles conj(Poles)],1);
                A(n,n,1:numel(a)-1)= -a(2:end);
            end

            % the second part is the Moving average (filters for connections)
            for n1 = 1:obj.NodeNum
                for n2 = 1:obj.NodeNum
                    if ~isempty(obj.Connections(n1,n2).FiltType) && (n1~=n2) 
                        lf = obj.Connections(n1,n2).LF; hf = obj.Connections(n1,n2).HF;
                        b1 = fir1(obj.Connections(n1,n2).FiltOrder, [lf hf]./(obj.SF/2),obj.Connections(n1,n2).FiltType);
                        A(n1,n2,1:numel(b1))= b1*obj.Connections(n1,n2).Gain/3;
                    end
                end
            end
            
            obj.ARMatrix = A;
        end
     
        %% Realization
        function [obj,TS] = Realization(obj,NS)
            alpha = 0.2;   
            order = size(obj.ARMatrix,3);
            if isempty(obj.ARMatrix)
               obj = GenerateARMatrix(obj);
            end
            TS = zeros(obj.NodeNum,order);
            for t = order+1:10000+NS
                TS_temp = zeros(obj.NodeNum,1);
                for ord = 1:order
                    TS_temp = TS_temp +obj.ARMatrix(:,:,ord)*TS(:,t-ord)+randn(1,1).*alpha;
                end
                TS(:,t) = TS_temp;
            end
            TS = TS(:,10001:end);% remove the first samples
            for n = 1:obj.NodeNum
                obj.Nodes(n).TS = TS(n,:);
            end
        end
        
        %% Stability check

    end
    
end

%% Other functions
function [r,f,Gain,F1,P1] = parameter_search(f,r0,P)
% set the amplitude of the poles in a way that gain(f1)/gain(f2)=P
% This is a recursive function

%initial values
f = sort(f); % sort frequencies from low to high
if ~exist('r0','var') || isempty(r0)
    r0 = .95;% = repmat(.95,size(f));
end
r = repmat(r0,size(f));

if numel(f)>1
    for p = numel(f)-1:-1:1
        %  set parameters for search
        P1= 0;
        e = 0.005;% confidence interval
        dis = 1-r(1);flag=0;
        while ((P1>P+e) || (P1<P-e))
            [XP,YP] = pol2cart([f(p:p+1) -f(p:p+1)],[r(p:p+1) r(p:p+1)]);% poles location
            [XT,YT] = pol2cart(f(p:p+1),ones(size(f(p:p+1)))); % the e^jw location
            Dists = squareform(pdist([XP' YP';XT' YT']));
            L = Dists(end-1:end,1:end-2);
            P1 = prod(L(2,:))./prod(L(1,:));
            % search the parameter space
            if flag == 1
                dis = dis/2;
            end
            if P1>P
                r(p) = (r(p)-dis);
            else
                flag=1;
                r(p) = r(p)+dis;
            end
        end
    end
end
F1 = 0:.02:pi;
for F= 1:numel(F1)
    [XP,YP] = pol2cart([f -f],[r r]);% poles location
    [XT,YT] = pol2cart(F1(F),1); % the e^jw location
    Dists = squareform(pdist([XP' YP';XT' YT']));
    L = Dists(end,1:2*numel(f));
    P1(F) = 1./prod(L(1,:));
end
Gain = max(P1);
end


