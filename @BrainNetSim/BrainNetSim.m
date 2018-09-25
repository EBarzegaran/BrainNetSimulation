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
        function obj = BrainNetSim(NodeNum,SF)
            obj.NodeNum = NodeNum;
            obj.Nodes = struct('Type',cell(NodeNum,1),'Freq',cell(NodeNum,1),'TS',cell(NodeNum,1));
            obj.Connections = struct('FiltType',cell(NodeNum),'FiltOrder',cell(NodeNum),'LF',cell(NodeNum),'HF',cell(NodeNum),'Gain',cell(NodeNum));
            if ~exist('SF','var') || isempty(SF)
                obj.SF = 100;
            end
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
            
            % 
            for n = 1: obj.NodeNum
                if numel(obj.Nodes(n).Freq)==1 % node with one oscillation
                    f1 = obj.Nodes(n).Freq/obj.SF; % first frequency
                    r1 = .95;
                    COS1 = cos(2*pi*f1);
                    A(n,n,1) = 2*r1*COS1;
                    A(n,n,2) = -r1^2;
                else % node with two oscillations
                    f1= obj.Nodes(n).Freq(1)/obj.SF; % first frequency
                    f2= obj.Nodes(n).Freq(2)/obj.SF; % second frequency
                    r1 = .85;%.87; % their relative amplitude
                    r2 = .9;%.96; 
                    COS1 = cos(2*pi*f1);
                    COS2 = cos(2*pi*f2);

                    A(n,n,1) = 2*r1*COS1 + 2*r2*COS2;
                    A(n,n,2) = -(r1^2 + r2^2 + 4*r1*r2*COS1*COS2);
                    A(n,n,3) = (2*r1*r2^2*COS1+2*r1^2*r2*COS2);
                    A(n,n,4) = -r1^2*r2^2;
                end
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
        
        %Stability check
        
        
    end
    
end
