function nodes_offset=offsetgeometry(nodes,offset,varargin)
%% Function to offset closed section
%
% Inputs:
% nodes: [N*2] matrix with 2D coordinates
% offset: vector with offset values for elements between nodes, if scalar then equal for all
%
% Outputs:
% nodes_offset: [N*2] matrix with coordinates with offset
%
%% Parse inputs

p=inputParser;
addParameter(p,'plot',false) % plot section
addParameter(p,'plottemp',false) % plot section
% addParameter(p,'offsetdir','in') % offset section outwards or inwards

parse(p,varargin{:});
do_plot=p.Results.plot;
do_plottemp=p.Results.plottemp;
% offsetdir=p.Results.offsetdir;

%%

% Homogeneous offset
if length(offset)==1
    offset=offset*ones(size(nodes,1),1);
end

if do_plot
    figure(); hold on; grid on; axis image;
    % plot(nodes(:,1),nodes(:,2),'ob');
end

nodes_offset=zeros(size(nodes));
n1_cell=cell(size(nodes,1),1);
n2_cell=cell(size(nodes,1),1);
N_cell=cell(size(nodes,1),1);
N_offset_cell=cell(size(nodes,1),1);

for k=1:length(nodes)
    
    N1=nodes(k,:);
    if k==length(nodes); ind=1; else; ind=k+1; end
    N2=nodes(ind,:);
    
    % n1 is vector along element
    n1=N2-N1;
    L=norm(n1);
    
    % Find tangent, and normal
    n1=n1/norm(n1);
    n3=[0 0 1];
    n2=cross(n3,[n1 0]); n2=n2(1:2)/norm(n2(1:2));
    
    n1_cell{k}=n1;
    n2_cell{k}=n2;
    
    % Offset element initially
    N1_offset_temp=N1+offset(k)*n2;
    N2_offset_temp=N2+offset(k)*n2;
    
    if do_plot
        h1(k)=plot([N1(1) N2(1)],[N1(2) N2(2)],'-ob','DisplayName','Initial');
    end
    
    if do_plottemp
        h2(k)=plot([N1_offset_temp(1) N2_offset_temp(1)],[N1_offset_temp(2) N2_offset_temp(2)],'-or','DisplayName','Offset temp');
    end
    
    N_cell{k}{1}=N1;
    N_cell{k}{2}=N2;
    
    N_offset_cell{k}{1}=N1_offset_temp;
    N_offset_cell{k}{2}=N2_offset_temp;
    
end

% Find intersection between new offset elements
for k=1:length(nodes)
    
    % Parametrized vectors:
    % v1_a=N1_a+t_a*n1_a
    % v1_b=N1_b+t_b*n1_b
    % Solve v1_a=v1_b
    
    N1_a=N_offset_cell{k}{1};
    n1_a=n1_cell{k};
    if k==length(nodes); ind=1; else; ind=k+1; end
    N1_b=N_offset_cell{ind}{1};
    n1_b=n1_cell{ind};
    
    [t_a,t_b,Intersect_a,Intersect_b]=LineIntersect(N1_a,N1_b,n1_a,n1_b);
    
    nodes_offset(k,:)=Intersect_a;
    
end

nodes_offset=[nodes_offset(end,:); nodes_offset(1:end-1,:)];

for k=1:length(nodes)
    
    N1=nodes_offset(k,:);
    if k==length(nodes); ind=1; else; ind=k+1; end
    N2=nodes_offset(ind,:);
    
    if do_plot
        h3(k)=plot([N1(1) N2(1)],[N1(2) N2(2)],'-xk','DisplayName','Offset');
    end
    
end

if do_plot & do_plottemp
    legend([h1(1) h2(1) h3(1)],'Location','NorthOutside');
    axistight(gca,[0.2 0.2],'x','y');
elseif do_plot
    legend([h1(1) h3(1)],'Location','NorthOutside');
    axistight(gca,[0.2 0.2],'x','y');
end


% axistight(gca,[0.2 0.2],'x','y');