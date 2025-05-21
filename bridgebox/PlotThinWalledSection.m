function varargout=plotthinwalledsection(nodes,elements,thickness,yc,zc,ys,zs,varargin)
%% Plot thin-walled section (box girder)
%
% nodes: matrix with rows [nodenumber, coord_y, coord_z]
% elements: matrix with rows [el_number, nodenumber1, nodenumber2]
% thickness: vector with element thickness corresponding to elements
% yc: center of area
% zc: center of area
% ys: center of shear
% zs: center of shear
%
%% Parse inputs

p=inputParser;
addParameter(p,'markersize',3)
addParameter(p,'nodecolor',[1 0 0])
addParameter(p,'elcolor',[0 0 1])
addParameter(p,'linewidth',0.75)
addParameter(p,'plottext',true,@islogical)
addParameter(p,'plotthickness',true,@islogical)
addParameter(p,'fontsize',6)
addParameter(p,'shift',false,@islogical)
addParameter(p,'offset',[],@isnumeric)
addParameter(p,'mirror',false,@islogical)
addParameter(p,'hax',[],@ishandle)
addParameter(p,'interpreter','latex',@ischar)
addParameter(p,'xlabel',[],@ischar)
addParameter(p,'ylabel',[],@ischar)

parse(p,varargin{:});
markersize=p.Results.markersize;
nodecolor=p.Results.nodecolor;
elcolor=p.Results.elcolor;
linewidth=p.Results.linewidth;
plottext=p.Results.plottext;
plotthickness=p.Results.plotthickness;
fontsize=p.Results.fontsize;
shift=p.Results.shift;
offset=p.Results.offset;
mirror=p.Results.mirror;
hax=p.Results.hax;
interpreter=p.Results.interpreter;
xl=p.Results.xlabel;
yl=p.Results.ylabel;

%%

% If no axes handle provided, create new figure
if isempty(hax)
    figure();
    ha=tight_subplot(1,1,[0.1],[0.05],0.05);
    sizefig();
else
    axes(hax);
end

hold on; grid on;
axis image;

% If elements is empty, assume a closed section with elements from node 1 to node 2
if isempty(elements)
    el_number=1:size(nodes,1);
    node_number_start=nodes(:,1);
    node_number_end=[ node_number_start(2:end) ; node_number_start(1)];
    elements=[ el_number.' node_number_start node_number_end];
end


% If no thickness provided,
if isempty(thickness)
    thickness=NaN*ones(size(elements,1),1);
end

% shift=true shifts the section to the center of area
if shift==true
    
    nodes_shift=nodes-[0 yc zc];
    yc_shift=yc-yc;
    zc_shift=yc-yc;
    ys_shift=ys-yc;
    zs_shift=zs-zc;
    
    nodes=nodes_shift;
    yc=yc_shift;
    zc=zc_shift;
    ys=ys_shift;
    zs=zs_shift;
    
end

% Shift whole section
if ~isempty(offset)
    
    nodes=nodes+[0 offset];
    yc=yc+offset(1);
    zc=zc+offset(2);
    
    ys=ys+offset(1);
    zs=zs+offset(2);
end

% Mirror section left-right
if mirror==true
    
    nodes1=nodes;
    nodes2=[nodes(:,1)+100e3 -nodes(:,2) nodes(:,3)];
    
    elements1=elements;
    elements2=elements+100e3;
    
    thickness1=thickness;
    thickness2=[thickness(:,1)+100e3];
    
    nodes=[nodes1 ; nodes2];
    elements=[elements1 ; elements2];
    thickness=[thickness1 ; thickness2];
    
end

% Plot nodes as markers
if markersize==0 | isnan(markersize)
else
h_nodes=plot(nodes(:,2),nodes(:,3),'o','markersize',markersize,'Color',nodecolor);
hideanno(h_nodes);
end


% Plot node labels
if plottext
    for k=1:length(nodes)
        text(nodes(k,2),nodes(k,3),['N' num2str(nodes(k,1))],'fontsize',fontsize);
    end
end


plot_matrix_x={};
plot_matrix_y={};

% Plot elements
for k=1:size(elements,1)
    
    % idxex of node 1 and node 2
    idx1=find(elements(k,2)==nodes(:,1));
    idx2=find(elements(k,3)==nodes(:,1));
    
    N1=nodes(idx1,2:3);
    N2=nodes(idx2,2:3);
    
    hp_el(k)=plot([N1(1) N2(1)],[N1(2) N2(2)],'-','Color',elcolor,'linewidth',linewidth);

    plot_matrix_x{end+1}=[N1(1) N2(1)];
    plot_matrix_x{end+1}=[NaN];

    plot_matrix_y{end+1}=[N1(2) N2(2)];
    plot_matrix_y{end+1}=[NaN];

    % Plot element labels
    if plottext
        
        element_text=['E' num2str(elements(k,1))];
        
        if plotthickness
            element_text=[element_text ', t=' num2str(thickness(k,1),'%0.1e')];
        end
        
        text(mean([N1(1) N2(1)]),mean([N1(2) N2(2)]),element_text,'fontsize',fontsize);
    end
    
end

delete(hp_el)

hp_el=plot(cell2mat(plot_matrix_x),cell2mat(plot_matrix_y),'-','Color',elcolor,'linewidth',linewidth);

% Plot center of area
if ~isempty(yc) & ~isempty(zc)
    plot(yc,zc,'xb');
    plot(yc,zc,'ob');
    text(yc,zc,'  CA','Color','b','fontsize',fontsize,'BackGroundColor','none');
end

% Plot shear center
if ~isempty(ys) & ~isempty(zs)
    plot(ys,zs,'xr');
    plot(ys,zs,'or');
    text(ys,zs,'CS  ','Color','r','fontsize',fontsize,'VerticalAlignment','Top','HorizontalAlignment','Right','BackGroundColor','none');
end

xlim([min(nodes(:,2))-1 max(nodes(:,2))+1]);
ylim([min(nodes(:,3))-0.5 max(nodes(:,3))+0.5]);
axistight(gca,[0.1 0.1],'x','y');

if isempty(xl); xl='$y$ [m]'; end
if isempty(yl); yl='$z$ [m]'; end


xlabel(xl,'Interpreter',interpreter);
ylabel(yl,'Interpreter',interpreter);

% if nargout>1
varargout{1}=hp_el;
% end
