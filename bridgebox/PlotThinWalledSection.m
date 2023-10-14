function PlotThinWalledSection(Nodes,Elements,Thickness,yc,zc,ys,zs,varargin)
%% Plot thin-walled section (box girder)
%
% Nodes: matrix with rows [nodenumber, coord_y, coord_z]
% Elements: matrix with rows [elnumber, nodenumber1, nodenumber2]
% Thickness: vector with element thickness corresponding to elements
% yc: center of area
% zc: center of area
% ys: center of shear
% zs: center of shear
%
%% Parse inputs

p=inputParser;
addParameter(p,'MarkerSize',3)
addParameter(p,'NodeColor',[1 0 0])
addParameter(p,'ElColor',[0 0 1])
addParameter(p,'LineWidth',0.75)
addParameter(p,'PlotText',true,@islogical)
addParameter(p,'PlotThickness',true,@islogical)
addParameter(p,'FontSize',6)
addParameter(p,'shift',false,@islogical)
addParameter(p,'offset',[],@isnumeric)
addParameter(p,'mirror',false,@islogical)
addParameter(p,'hax',[],@ishandle)

parse(p,varargin{:});
MarkerSize=p.Results.MarkerSize;
NodeColor=p.Results.NodeColor;
ElColor=p.Results.ElColor;
LineWidth=p.Results.LineWidth;
PlotText=p.Results.PlotText;
PlotThickness=p.Results.PlotThickness;
FontSize=p.Results.FontSize;
shift=p.Results.shift;
offset=p.Results.offset;
mirror=p.Results.mirror;
hax=p.Results.hax;

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
if isempty(Elements)
    ElNumber=1:size(Nodes,1);
    NodeNumberStart=Nodes(:,1);
    NodeNumberEnd=[ NodeNumberStart(2:end) ; NodeNumberStart(1)];
    Elements=[ ElNumber.' NodeNumberStart NodeNumberEnd];
end


% If no thickness provided,
if isempty(Thickness)
    Thickness=NaN*ones(size(Elements,1),1);
end

% shift=true shifts the section to the center of area
if shift==true
    
    Nodes_shift=Nodes-[0 yc zc];
    yc_shift=yc-yc;
    zc_shift=yc-yc;
    ys_shift=ys-yc;
    zs_shift=zs-zc;
    
    Nodes=Nodes_shift;
    yc=yc_shift;
    zc=zc_shift;
    ys=ys_shift;
    zs=zs_shift;
    
end

% Shift whole section
if ~isempty(offset)
    
    Nodes=Nodes+[0 offset];
    yc=yc+offset(1);
    zc=zc+offset(2);
    
    ys=ys+offset(1);
    zs=zs+offset(2);
end

% Mirror section left-right
if mirror==true
    
    Nodes1=Nodes;
    Nodes2=[Nodes(:,1)+100e3 -Nodes(:,2) Nodes(:,3)];
    
    Elements1=Elements;
    Elements2=Elements+100e3;
    
    Thickness1=Thickness;
    Thickness2=[Thickness(:,1)+100e3];
    
    Nodes=[Nodes1 ; Nodes2];
    Elements=[Elements1 ; Elements2];
    Thickness=[Thickness1 ; Thickness2];
    
end

% Plot nodes are markers
plot(Nodes(:,2),Nodes(:,3),'o','MarkerSize',MarkerSize,'Color',NodeColor);

% Plot node labels
if PlotText
    for k=1:length(Nodes)
        text(Nodes(k,2),Nodes(k,3),['N' num2str(Nodes(k,1))],'FontSize',FontSize);
    end
end

% Plot elements
for k=1:size(Elements,1)
    
    % Index of node 1 and node 2
    ind1=find(Elements(k,2)==Nodes(:,1));
    ind2=find(Elements(k,3)==Nodes(:,1));
    
    N1=Nodes(ind1,2:3);
    N2=Nodes(ind2,2:3);
    
    plot([N1(1) N2(1)],[N1(2) N2(2)],'-','Color',ElColor);
    
    % Plot element labels
    if PlotText
        
        ElementText=['E' num2str(Elements(k,1))];
        
        if PlotThickness
            ElementText=[ElementText ', t=' num2str(Thickness(k,1),'%0.1e')];
        end
        
        text(mean([N1(1) N2(1)]),mean([N1(2) N2(2)]),ElementText,'FontSize',FontSize);
    end
    
end

% Plot center of area
if ~isempty(yc) & ~isempty(zc)
    plot(yc,zc,'xb');
    plot(yc,zc,'ob');
    text(yc,zc,'  CA','Color','b','FontSize',FontSize,'BackGroundColor','none');
end

% Plot shear center
if ~isempty(ys) & ~isempty(zs)
    plot(ys,zs,'xr');
    plot(ys,zs,'or');
    text(ys,zs,'CS  ','Color','r','FontSize',FontSize,'VerticalAlignment','Top','HorizontalAlignment','Right','BackGroundColor','none');
end

xlim([min(Nodes(:,2))-1 max(Nodes(:,2))+1]);
ylim([min(Nodes(:,3))-0.5 max(Nodes(:,3))+0.5]);
axistight(gca,[0.1 0.1],'x','y');

xlabel('y [m]');
ylabel('z [m]');

