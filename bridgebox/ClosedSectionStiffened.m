function [Nodes,Elements,Thickness]=ClosedSectionStiffened(NodesBox,ThicknessBox,StiffenerType,StiffenerCenterDist,StiffenerGeo,StiffenerThickness,varargin)
%% Design (closed) bridge section with stiffeners
%
% Inputs:
% NodesBox: [N,2] matrix with rows [Coord_y, Coord_z] for box shape
% ThicknessBox: [1,N] vector with outer skin thickness for each segment
% StiffenerType: [1,N] vector with stiffener types for each segment
% StiffenerCenterDist: [1,N] vector with stiffener cc distance for each segment
% StiffenerGeo: [1,M] cell with geometry for each type of stiffener, e.g. trapezoidal, rectangular, or knife
% StiffenerThickness: [1,M] vector with skin thickness for each type of stiffener
%
% Outputs:
% Nodes: [*,3] matrix with rows [NodeNumber Coord_y Coord_z]
% Elements: [*,3] matrix with rows [ElementNumber NodeNumber1 NodeNumber2]
% Thickness: [*,1] matrix with rows [Thickness] corresponding to the elements
%
% Variable inputs:
% If the section is already designed, use options DistanceEdgeStart and N_stiff to get the correct design
% For a new design, use option DistanceEdgeMin to place stiffeners automatically
%
%% Parse inputs

p=inputParser;
addParameter(p,'Plot',true)
addParameter(p,'PlotStiffener',false)
addParameter(p,'PlotText',true)
addParameter(p,'DistanceEdgeMin',0.1*ones(size(StiffenerCenterDist)))
addParameter(p,'DistanceEdgeStart',[])
addParameter(p,'N_stiff',[])
addParameter(p,'tol',10e-3)

parse(p,varargin{:});
DoPlot=p.Results.Plot;
DoPlotStiffener=p.Results.PlotStiffener;
PlotText=p.Results.PlotText;
DistanceEdgeMin=p.Results.DistanceEdgeMin;
DistanceEdgeStart=p.Results.DistanceEdgeStart;
N_stiff=p.Results.N_stiff;
tol=p.Results.tol;

%% Check inputs

% Ensure correct dimension
if size(NodesBox,2)>size(NodesBox,1)
    NodesBox=NodesBox.';
end

% Number of elements in the box (=number of nodes since closed section)
N_el=size(NodesBox,1);

if length(ThicknessBox)~=N_el
    error('ThicknessBox must have N_el elements');
end

if length(StiffenerType)~=N_el
    error('StiffenerType must have N_el elements');
end

if length(StiffenerCenterDist)~=N_el
    error('StiffenerCenterDist must have N_el elements');
end

if length(DistanceEdgeMin)~=N_el
    error('DistanceEdgeMin must have N_el elements');
end

if length(StiffenerGeo)~=length(StiffenerThickness)
    error('StiffenerGeo must be same size as StiffenerThickness');
end

if ~isempty(DistanceEdgeStart)
    if isempty(N_stiff)
        error('Number of stiffeners (N_stiff) must be non-empty if DistanceEdgeStart is given');
    end
end

% If the signed area is negative, the nodes are not in CCW order
[A_signed]=PolygonMoments(NodesBox);
if A_signed<0
    NodesBox=flip(NodesBox,1);
end

%%

% Find coordinates of start/end node for each box segment
for k=1:size(NodesBox,1)
    
    if k==size(NodesBox,1); k_next=1;
    else; k_next=k+1;
    end
    
    N1{k}=NodesBox(k,:);
    N2{k}=NodesBox(k_next,:);
    L_el(k)=norm(N2{k}-N1{k});
    
end

% Check that stiffener start in (0,0)
for k=1:length(StiffenerGeo)

    if any(isnan(StiffenerGeo{k})) | isempty(StiffenerGeo{k})
        continue
    end

    if StiffenerGeo{k}(1,1)~=0 | StiffenerGeo{k}(1,2)~=0
        k
        StiffenerGeo{k}
        error('First coordinate of stiffener must be (0,0)');
	end
end

% Find length of base of each stiffener type
for k=1:length(StiffenerGeo)

    % If no stiffener, continue
    if isnan(StiffenerType(k))
        StiffenerBase(k)=NaN;
        continue
    end
    
    % If z-coordinate of last point is not zero, the stiffener is open (e.g. knife)
    if StiffenerGeo{k}(end,2)~=0
        StiffenerBase(k)=0;
    else
        StiffenerBase(k)=norm(StiffenerGeo{k}(1,:)-StiffenerGeo{k}(end,:));
    end
end

% Save nodes of box corners
Nodes_corner=[ [1:length(NodesBox)].' NodesBox];

for k=1:size(NodesBox,1)
    
    % Initiate
    Nodes_outer{k}=[];
    El_outer{k}=[];
    Thickness_outer{k}=[];
    
    Nodes_stiff{k}=[];
    El_stiff{k}=[];
    Thickness_stiff{k}=[];
    
    % Base number
    NodeNumberOuter=1e2*k;
    ElementNumberOuter=1e2*k;
    
    NodeNumberStiff=5e3+1e2*k;
    ElementNumberStiff=5e3+1e2*k;
    
    if k==size(NodesBox,1); k_next=1;
    else; k_next=k+1;
    end
    
    % If no stiffener, then just add element from corner to corner
    if isnan(StiffenerType(k))
        ElementNumberOuter=ElementNumberOuter+1;
        El_outer{k}(end+1,:)=[ElementNumberOuter k k_next];
        Thickness_outer{k}(end+1,:)=[ElementNumberOuter ThicknessBox(k)];
        continue
    end
    
    % Base distance of stiffener used for this element
    DistanceBase=StiffenerBase(StiffenerType(k));
    
    % Determine number of stiffeners, if DistanceEdgeStart not provided then design new, else add N_stiff
    if isempty(DistanceEdgeStart)
        % Add stiffeners until DistanceEnd is smaller than DistanceEdgeMin
        for n_stiff=0:100
            DistanceEdge=(L_el(k)-(n_stiff-1)*StiffenerCenterDist(k)-DistanceBase/2*2)/2;
            if DistanceEdge<DistanceEdgeMin(k)
                n_stiff=n_stiff-1;
                DistanceEdge=(L_el(k)-(n_stiff-1)*StiffenerCenterDist(k))/2;
                
                if n_stiff==0
                    warning(['No stiffener added to segment ' num2str(k) ', space too small']);
                end
                
                break
            end
        end
        
    else
        DistanceEdge=DistanceEdgeStart(k);
        n_stiff=N_stiff(k);
    end
    
    % Vector along box element
    n1=N2{k}-N1{k};
    
    % Rotation angle of element relative to x-axis and rotation matrix
    angle_rad=atan2(n1(2),n1(1));
    T=[cos(angle_rad) -sin(angle_rad) ; sin(angle_rad) cos(angle_rad)];
    
    % Add stiffeners
    for j=1:n_stiff
        
        if j==1; NodeNumberOuterPrev=k;
        else; NodeNumberOuterPrev=NodeNumberOuter;
        end
        
        % Starting node stiffener
        NodeNumberOuter=NodeNumberOuter+1;
        t_a=(DistanceEdge+(j-1)*StiffenerCenterDist(k)-DistanceBase/2)/L_el(k);
        Nodes_outer{k}(end+1,:)=[NodeNumberOuter N1{k}+n1*t_a];
        
        % End node stiffener
        NodeNumberOuter=NodeNumberOuter+1;
        t_b=(DistanceEdge+(j-1)*StiffenerCenterDist(k)+DistanceBase/2)/L_el(k);
        Nodes_outer{k}(end+1,:)=[NodeNumberOuter N1{k}+n1*t_b];
        
        if t_a>1
            k
            j
            warning('Stiffener outside element');
        end
        
        if t_b>1
            k
            j
            warning('Stiffener outside element');
        end
        
        % Element before stiffener
        ElementNumberOuter=ElementNumberOuter+1;
        El_outer{k}(end+1,:)=[ElementNumberOuter NodeNumberOuterPrev NodeNumberOuter-1];
        Thickness_outer{k}(end+1,:)=[ElementNumberOuter ThicknessBox(k)];
        
        % Element base of stiffener
        ElementNumberOuter=ElementNumberOuter+1;
        El_outer{k}(end+1,:)=[ElementNumberOuter NodeNumberOuter-1 NodeNumberOuter];
        Thickness_outer{k}(end+1,:)=[ElementNumberOuter ThicknessBox(k)];
        
        % Nodes and elements or stiffener body
        N_stiff_nodes=size(StiffenerGeo{StiffenerType(k)},1);
        for h=1:(N_stiff_nodes)
            
            NodeNumberStiff=NodeNumberStiff+1;
            Nodes_stiff{k}(end+1,:)=[NodeNumberStiff (T*StiffenerGeo{StiffenerType(k)}(h,:).').'+N1{k}+n1*t_a];
            
            % For second node and on, add element
            if h>1
                ElementNumberStiff=ElementNumberStiff+1;
                El_stiff{k}(end+1,:)=[ElementNumberStiff NodeNumberStiff-1 NodeNumberStiff];
                Thickness_stiff{k}(end+1,:)=[ElementNumberStiff StiffenerThickness(StiffenerType(k))];
            end
            
        end
        
    end
    
    % Element after stiffeners
    ElementNumberOuter=ElementNumberOuter+1;
    El_outer{k}(end+1,:)=[ElementNumberOuter NodeNumberOuter k_next];
    Thickness_outer{k}(end+1,:)=[ElementNumberOuter ThicknessBox(k)];
    
end

%% Merge

Nodes=stackVertical(Nodes_corner,Nodes_outer{:},Nodes_stiff{:});
Elements=stackVertical(El_outer{:},El_stiff{:});
Thickness=stackVertical(Thickness_outer{:},Thickness_stiff{:});

Thickness=Thickness(:,2);

%% Delete nodes within tolerance

dy=Nodes(:,2)-Nodes(:,2).';
dz=Nodes(:,3)-Nodes(:,3).';

Z=sqrt(dy.^2+dz.^2);

ind_delete=[];
for ind1=1:size(Z,1)
    for ind2=1:(ind1-1)
        
        % Merge nodes within distance tol, box outer nodes always master
        if Z(ind1,ind2)<tol
            
            Node_master=Nodes(ind2,1);
            Node_slave=Nodes(ind1,1);
            
            Nodes(ind1,:)=NaN;
            TempNodes=Elements(:,2:3);
            TempNodes(TempNodes==Node_slave)=Node_master;
            Elements(:,2:3)=TempNodes;
            
            ind_delete=[ind_delete ind1];
        end
        
    end
end
Nodes(ind_delete,:)=[];

% Delete elements from node A to node A (zero length)
ind_delete=[];
for k=1:length(Elements)
    if Elements(k,2)==Elements(k,3)
        ind_delete=[ind_delete k];
    end
end
Elements(ind_delete,:)=[];
Thickness(ind_delete,:)=[];

%% Plot

if DoPlot
    PlotThinWalledSection(Nodes,Elements,Thickness,[],[],[],[],'PlotText',PlotText,'FontSize',8);
end

if DoPlotStiffener
    
    figure(); 
    
    ha=tight_subplot(1,length(StiffenerGeo),0.1,[0.1 0.1],[0.1 0.05]);
    
    for k=1:length(StiffenerGeo)
        StiffNodes=[ [1:size(StiffenerGeo{k},1)].' StiffenerGeo{k}];
        StiffElements=[  [1:(size(StiffNodes,1)-1)].' [1:(size(StiffNodes,1)-1)].'  [2:size(StiffNodes,1)].' ];
        StiffThickness=[StiffenerThickness(k)];
        PlotThinWalledSection(StiffNodes,StiffElements,Thickness,[],[],[],[],'PlotText',PlotText,'FontSize',8,'hax',ha(k));
        
    end
end

end

