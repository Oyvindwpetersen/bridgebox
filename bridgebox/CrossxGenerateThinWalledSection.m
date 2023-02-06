function CrossxGenerateThinWalledSection(Nodes,Elements,Thickness,name_section,name_file,plot_section,varargin)

%% Generate CrossX input file for thinwalled section
%
% Inputs:
% Nodes: (n_node*3) matrix with rows: node number, x coord, y coord
% Elements: (n_el*3) matrix with rows: el number, start node number, end node number
% thickness_vector: vector with thickness of elements
% name_section: name of section
% name_file: name of txt file
% plot_section: logical argument (true/false) to plot
%
% Outputs:
% saved txt file that can be imported in crossx
%
%% Parse inputs

p=inputParser;
addParameter(p,'unit','mm') % Input units
addParameter(p,'MarkerSize',3)
addParameter(p,'NodeColor',[1 0 0])
addParameter(p,'ElColor',[0 0 1])
addParameter(p,'PlotText',true)
addParameter(p,'FontSize',6)

parse(p,varargin{:});
unit=p.Results.unit;
MarkerSize=p.Results.MarkerSize;
NodeColor=p.Results.NodeColor;
ElColor=p.Results.ElColor;
PlotText=p.Results.PlotText;
FontSize=p.Results.FontSize;

%%

if strcmpi(unit,'m')
    scale_factor=1000;
else
    scale_factor=1;
end

Nodes(:,2:3)=scale_factor*Nodes(:,2:3);
Thickness=scale_factor*Thickness;

n_el=size(Elements,1);
n_node=size(Nodes,1);

if length(Thickness)==1
    Thickness=Thickness*ones(n_el,1);
end

if length(Thickness)~=n_el
    error(['thickness_vector wrong length' '(is ' num2str(length(Thickness)) ', should be ' num2str(n_el) ')' ]);
end

n_matrix=Nodes(:,1);
xy_matrix=Nodes(:,2:3);

% Format matrix for CrossX, [node1_x node1_y node2_x node2_y ];
for k=1:n_el
    
    n1=find(Elements(k,2)==n_matrix);
    n2=find(Elements(k,3)==n_matrix);
    
    if isempty(n1)
        error(['First node of element ' num2str(Elements(k,1)) 'not found' ]);
    elseif isempty(n2)
        error(['Second node of element ' num2str(Elements(k,1)) 'not found' ]);
    end
    
    matrix_out(k,:)=[ [xy_matrix(n1,1) xy_matrix(n1,2) xy_matrix(n2,1) xy_matrix(n2,2)] Thickness(k)];
    

    L_el=sqrt( (matrix_out(k,1)-matrix_out(k,3))^2 + (matrix_out(k,2)-matrix_out(k,4))^2);

    
    if Thickness(k)/L_el > 0.5

        Thickness(k)
        L_el
        k
        warning(['Thickness of element more than 50% of element length']);
    end


end

%%

% Print txt file for CrossX
if ~isempty(name_file)
    
    formatSpecString = '%s \n';
    formatSpecNumEng = '%0.5e \n';
    
    clear inputData
    inputData{1,1}='§Section';
    inputData{2,1}=['# ' name_section];
    
    for k=1:n_el
        
        isnegative=matrix_out(k,:)<0;
        
        inputData{end+1,1}=[ repmat(' ',1,3-isnegative(1)) ... ...
            num2str(matrix_out(k,1),formatSpecNumEng) repmat(' ',1,7-isnegative(2)) ...
            num2str(matrix_out(k,2),formatSpecNumEng) repmat(' ',1,7-isnegative(3)) ...
            num2str(matrix_out(k,3),formatSpecNumEng) repmat(' ',1,7-isnegative(4)) ...
            num2str(matrix_out(k,4),formatSpecNumEng) repmat(' ',1,7-isnegative(5)) ...
            num2str(matrix_out(k,5),formatSpecNumEng)  ...
            ];
    end
    inputData{end+1,1}='§End';
    
    fid_input = fopen([name_file '.txt'],'wt');
    for k=1:size(inputData,1)
        fprintf(fid_input,formatSpecString,inputData{k,:});
    end
    fclose(fid_input);
    
end

%%

if plot_section

    PlotThinWalledSection(Nodes,Elements,Thickness,[],[],[],[]);

end

