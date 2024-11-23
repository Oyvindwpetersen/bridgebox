function [nodes,elements,thickness]=closedsectionstiffened2(nodes_box,t_box,stiff_geo,stiff_t,varargin)
%% Design (closed) bridge section with stiffeners
%
% Inputs:
% nodes_box: [N,2] matrix with rows [Coord_y, Coord_z] for box shape
% t_box: [1,N] vector with outer skin thickness for each segment
% stiff_type: [1,N] vector with stiffener types for each segment
% stiff_cc: [1,N] vector with stiffener cc distance for each segment
% stiff_geo: [1,M] cell with geometry for each type of stiffener, e.g. trapezoidal, rectangular, or knife
% stiff_t: [1,M] vector with skin thickness for each type of stiffener
%
% Outputs:
% nodes: [*,3] matrix with rows [NodeNumber Coord_y Coord_z]
% elements: [*,3] matrix with rows [ElementNumber NodeNumber1 NodeNumber2]
% thickness: [*,1] matrix with rows [thickness] corresponding to the elements
%
% Variable inputs:
% If the section is already designed, use options distedgestart and N_stiff to get the correct design
% For a new design, use option distedgemin to place stiffeners automatically
%
%% Parse inputs

p=inputParser;
addParameter(p,'Plot',true,@islogical)
addParameter(p,'PlotStiffener',false,@islogical)
addParameter(p,'PlotText',true,@islogical)
% addParameter(p,'distedgemin',0.1*ones(size(stiff_cc)))
% addParameter(p,'distedgestart',[])
% addParameter(p,'modestiff','auto',@islogical)
% addParameter(p,'diststiff',{})
% addParameter(p,'typestiff',{})
% addParameter(p,'N_stiff',[])
addParameter(p,'edge_dist',[],@isnumeric)
addParameter(p,'stiff_dist_cc',{},@iscell)
addParameter(p,'stiff_type',[],@iscell)
addParameter(p,'tol',10e-3)


% edge_dist(k)=...;
% stiff_dist_cc{k}(j)=...;
% stiff_dist_base{k}(j)=...;
% stiff_type{k}(j)=...



parse(p,varargin{:});
DoPlot=p.Results.Plot;
DoPlotStiffener=p.Results.PlotStiffener;
PlotText=p.Results.PlotText;
% stiff=p.Results.stiff;
% distedgemin=p.Results.distedgemin;
% distedgestart=p.Results.distedgestart;
% N_stiff=p.Results.N_stiff;

edge_dist=p.Results.edge_dist;
stiff_dist_cc=p.Results.stiff_dist_cc;
stiff_type=p.Results.stiff_type;

tol=p.Results.tol;

%% Check inputs

% Ensure correct dimension
if size(nodes_box,2)>size(nodes_box,1)
    nodes_box=nodes_box.';
end

% Number of elements in the box (=number of nodes since closed section)
N_el=size(nodes_box,1);

if length(t_box)~=N_el
    error('t_box must have N_el elements');
end

if length(stiff_type)~=N_el
    error('stiff_type must have N_el elements');
end

% if length(stiff_cc)~=N_el
%     error('stiff_cc must have N_el elements');
% end

% if length(distedgemin)~=N_el
%     error('distedgemin must have N_el elements');
% end
% 
% if length(stiff_geo)~=length(stiff_t)
%     error('stiff_geo must be same size as stiff_t');
% end

% if ~isempty(distedgestart)
%     if isempty(N_stiff)
%         error('Number of stiffeners (N_stiff) must be non-empty if distedgestart is given');
%     end
% end

% If the signed area is negative, the nodes are not in CCW order
[A_signed]=PolygonMoments(nodes_box);
if A_signed<0
    nodes_box=flip(nodes_box,1);
end

%%

% Find coordinates of start/end node for each box segment
for k=1:size(nodes_box,1)
    
    if k==size(nodes_box,1); k_next=1;
    else; k_next=k+1;
    end
    
    N1{k}=nodes_box(k,:);
    N2{k}=nodes_box(k_next,:);
    L_el(k)=norm(N2{k}-N1{k});
    
end

% Check that stiffener start in (0,0)
for k=1:length(stiff_geo)

    if any(isnan(stiff_geo{k})) | isempty(stiff_geo{k})
        continue
    end

    if stiff_geo{k}(1,1)~=0 | stiff_geo{k}(1,2)~=0
        k
        stiff_geo{k}
        error('First coordinate of stiffener must be (0,0)');
	end
end

% Find length of base of each stiffener type
for k=1:length(stiff_geo)

    % If no stiffener, continue
    % if isempty(stiff_type(k))
    %     stiff_base(k)=NaN;
    %     continue
    % end
    
    % If z-coordinate of last point is not zero, the stiffener is open (e.g. knife)
    if stiff_geo{k}(end,2)~=0
        stiff_base(k)=0;
    else
        stiff_base(k)=norm(stiff_geo{k}(1,:)-stiff_geo{k}(end,:));
    end
end

% if strcmpi(modestiff,'manual')

% for k=1:size(nodes_box)

        % Determine number of stiffeners, if distedgestart not provided then design new, else add N_stiff
    % if isempty(distedgestart)
    %     % Add stiffeners until DistanceEnd is smaller than distedgemin
    %     for n_stiff=0:100
    % 
    % 
    % 
    % 
    % 
    %         dist_edge=(L_el(k)-(n_stiff-1)*stiff_cc(k)-dist_base/2*2)/2;
    %         if dist_edge<distedgemin(k)
    %             n_stiff=n_stiff-1;
    %             dist_edge=(L_el(k)-(n_stiff-1)*stiff_cc(k))/2;
    % 
    %             if n_stiff==0
    %                 warning(['No stiffener added to segment ' num2str(k) ', space too small']);
    %             end
    % 
    %             break
    %         end
    %     end
    % 
    % else
        % dist_edge=distedgestart(k);
        % n_stiff=N_stiff(k);


    % end

% end

% 
% edge_dist(k)=...;
% stiff_dist_cc{k}(j)=...;
% stiff_dist_base{k}(j)=...;
% stiff_type{k}(j)=...

for k=1:size(nodes_box)

    for j=1:length(stiff_type{k})
        idx_type=stiff_type{k}(j);
        stiff_dist_base{k}(j)=stiff_base(idx_type);
    end

end

%%
% Save nodes of box corners
nodes_corner=[ [1:length(nodes_box)].' nodes_box];

for k=1:size(nodes_box,1)
    
    % Initiate
    nodes_outer{k}=[];
    el_outer{k}=[];
    t_outer{k}=[];
    
    nodes_stiff_tmp{k}=[];
    el_stiff{k}=[];
    thickness_stiff{k}=[];
    
    % Base number
    node_number_outer=1e2*k;
    el_number_outer=1e2*k;
    
    node_number_stiff=5e3+1e2*k;
    el_number_stiff=5e3+1e2*k;
    
    if k==size(nodes_box,1); k_next=1;
    else; k_next=k+1;
    end
    
    % If no stiffener, then just add element from corner to corner
    if isempty(stiff_type{k})
        el_number_outer=el_number_outer+1;
        el_outer{k}(end+1,:)=[el_number_outer k k_next];
        t_outer{k}(end+1,:)=[el_number_outer t_box(k)];
        continue
    end
    
    % Vector along box element
    n1=N2{k}-N1{k};
    
    % Rotation angle of element relative to x-axis and rotation matrix
    angle_rad=atan2(n1(2),n1(1));
    T=[cos(angle_rad) -sin(angle_rad) ; sin(angle_rad) cos(angle_rad)];
    
    % Add stiffeners
    n_stiff=length(stiff_type{k});
    for j=1:n_stiff
        
        if j==1; node_number_outer_prev=k;
        else; node_number_outer_prev=node_number_outer;
        end

        % Base distance of stiffener used for this element
        % dist_base=stiff_base(stiff_type(k));

        % Starting node stiffener
        node_number_outer=node_number_outer+1;
        t_a=(edge_dist(k)...
             +stiff_dist_base{k}(1)/2*0 ...
             +sum(stiff_dist_cc{k}(1:j-1)) ...
             -stiff_dist_base{k}(j)/2) ...
             /L_el(k);
        % t_a=(dist_edge+(j-1)*stiff_cc(k)-dist_base/2)/;
        nodes_outer{k}(end+1,:)=[node_number_outer N1{k}+n1*t_a];
        
        % End node stiffener
        node_number_outer=node_number_outer+1;
        % t_b=(dist_edge+(j-1)*stiff_cc(k)+dist_base/2)/L_el(k);
        t_b=t_a+(stiff_dist_base{k}(j))/L_el(k);
        nodes_outer{k}(end+1,:)=[node_number_outer N1{k}+n1*t_b];
        
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
        el_number_outer=el_number_outer+1;
        el_outer{k}(end+1,:)=[el_number_outer node_number_outer_prev node_number_outer-1];
        t_outer{k}(end+1,:)=[el_number_outer t_box(k)];
        
        % Element base of stiffener
        el_number_outer=el_number_outer+1;
        el_outer{k}(end+1,:)=[el_number_outer node_number_outer-1 node_number_outer];
        t_outer{k}(end+1,:)=[el_number_outer t_box(k)];
        
        % nodes and elements or stiffener body
        N_stiff_nodes=size(stiff_geo{stiff_type{k}(j)},1);
        for h=1:(N_stiff_nodes)
            
            node_number_stiff=node_number_stiff+1;
            nodes_stiff_tmp{k}(end+1,:)=[node_number_stiff (T*stiff_geo{stiff_type{k}(j)}(h,:).').'+N1{k}+n1*t_a];
            
            % For second node and on, add element
            if h>1
                el_number_stiff=el_number_stiff+1;
                el_stiff{k}(end+1,:)=[el_number_stiff node_number_stiff-1 node_number_stiff];
                thickness_stiff{k}(end+1,:)=[el_number_stiff stiff_t(stiff_type{k}(j))];
            end
            
        end
        
    end
    
    % Element after stiffeners
    el_number_outer=el_number_outer+1;
    el_outer{k}(end+1,:)=[el_number_outer node_number_outer k_next];
    t_outer{k}(end+1,:)=[el_number_outer t_box(k)];
    
end

%% Merge

nodes=stackVertical(nodes_corner,nodes_outer{:},nodes_stiff_tmp{:});
elements=stackVertical(el_outer{:},el_stiff{:});
thickness=stackVertical(t_outer{:},thickness_stiff{:});

thickness=thickness(:,2);

%% Delete nodes within tolerance

dy=nodes(:,2)-nodes(:,2).';
dz=nodes(:,3)-nodes(:,3).';

Z=sqrt(dy.^2+dz.^2);

ind_delete=[];
for ind1=1:size(Z,1)
    for ind2=1:(ind1-1)
        
        % Merge nodes within distance tol, box outer nodes always master
        if Z(ind1,ind2)<tol
            
            node_master=nodes(ind2,1);
            node_slave=nodes(ind1,1);
            
            nodes(ind1,:)=NaN;
            Tempnodes=elements(:,2:3);
            Tempnodes(Tempnodes==node_slave)=node_master;
            elements(:,2:3)=Tempnodes;
            
            ind_delete=[ind_delete ind1];
        end
        
    end
end
nodes(ind_delete,:)=[];

% Delete elements from node A to node A (zero length)
ind_delete=[];
for k=1:length(elements)
    if elements(k,2)==elements(k,3)
        ind_delete=[ind_delete k];
    end
end
elements(ind_delete,:)=[];
thickness(ind_delete,:)=[];

%% Plot

if DoPlot
    plotthinwalledsection(nodes,elements,thickness,[],[],[],[],'PlotText',PlotText,'FontSize',8);
end

if DoPlotStiffener
    
    figure(); 
    
    ha=tight_subplot(1,length(stiff_geo),0.1,[0.1 0.1],[0.1 0.05]);
    
    for k=1:length(stiff_geo)
        nodes_stiff_tmp=[ [1:size(stiff_geo{k},1)].' stiff_geo{k}];
        elements_stiff_tmp=[  [1:(size(nodes_stiff_tmp,1)-1)].' [1:(size(nodes_stiff_tmp,1)-1)].'  [2:size(nodes_stiff_tmp,1)].' ];
        t_stiff_tmp=[stiff_t(k)]*ones(size(elements_stiff_tmp,1),1);
        plotthinwalledsection(nodes_stiff_tmp,elements_stiff_tmp,t_stiff_tmp,[],[],[],[],'PlotText',PlotText,'FontSize',8,'hax',ha(k));
        
    end
end

end

