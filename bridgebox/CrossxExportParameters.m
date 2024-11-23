function [par_val,parameter_names,par_struct]=crossxexportparameters(name_file)

%%

if exist(name_file)~=2
    error(['File ' name_files ' does not exist']);
end

%%
rehash

[fid_input,errormsg]=fopen(name_file, 'r');

if ~isempty(errormsg)
    disp('Error in reading input file:');
    disp(errormsg);
	return
end

data=textscan(fid_input, '%s', 'delimiter', '\n', 'whitespace', '');
data=data{1}; %read inputfile
fclose(fid_input);

parameter_names_long={...
    'Area (true)'
    'y'' of E-weighted c.o.g.'
    'z'' of E-weighted c.o.g.'
    'y'' of shear-center'
    'z'' of shear-center'
    'I_y'' (y'' at c.o.g.)'
    'I_z'' (z'' at c.o.g.)'
    'I_y''z'' (at c.o.g.)'
    'I_t'
    'I_omega (warping/Cw)'
    'I_r (polar)'
};

parameter_names={...
    'A'
    'y_prime_cog'
    'z_prime_cog'
    'y_prime_cs'
    'z_prime_cs'
    'I_y_prime_cog'
    'I_z_prime_cog'
    'I_yz_prime_cog'
    'I_t'
    'I_omega'
    'I_r'
};

% Scale from mm to m

scalefactor(1)=(1e3)^2; % Area
scalefactor(2)=(1e3)^1; % Distance
scalefactor(3)=(1e3)^1; % Distance
scalefactor(4)=(1e3)^1; % Distance
scalefactor(5)=(1e3)^1; % Distance
scalefactor(6)=(1e3)^4; % Second moment of area
scalefactor(7)=(1e3)^4; % Second moment of areapar_mm
scalefactor(8)=(1e3)^4; % Second moment of area
scalefactor(9)=(1e3)^4; % Torsional constant
scalefactor(10)=(1e3)^6; % Cw constant
scalefactor(11)=(1e3)^4; % Ir constant

par_val=[];
for k=1:length(parameter_names_long)

    line_idx=strmatch(lower(parameter_names_long{k}),lower(data));
    
    if isnan(line_idx); par_val(k,1)=NaN; continue; end
    
    char_idx=strfind(data{line_idx},':');
    % par_val(k,1)=inpReadNumber(data,line_idx,'%f','skip',char_idx);

    skip=char_idx;
    line=data{line_idx};

    par_val(k,1)=sscanf(line(skip+1:end),'%f')';

    par_val(k,1)=par_val(k,1)/scalefactor(k);
    
end

par_struct=struct();

for k=1:length(parameter_names_long)
	par_struct=setfield(par_struct,parameter_names{k},par_val(k));
end





