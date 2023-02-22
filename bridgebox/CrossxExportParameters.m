function [M,ParNamesRed,M_out]=CrossxExportParameters(name_file)

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

ParNames={...
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

ParNamesRed={...
    'Area'
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

ScaleFactor(1)=(1e3)^2; % Area
ScaleFactor(2)=(1e3)^1; % Distance
ScaleFactor(3)=(1e3)^1; % Distance
ScaleFactor(4)=(1e3)^1; % Distance
ScaleFactor(5)=(1e3)^1; % Distance
ScaleFactor(6)=(1e3)^4; % Second moment of area
ScaleFactor(7)=(1e3)^4; % Second moment of area
ScaleFactor(8)=(1e3)^4; % Second moment of area
ScaleFactor(9)=(1e3)^4; % Torsional constant
ScaleFactor(10)=(1e3)^6; % Cw constant
ScaleFactor(11)=(1e3)^4; % Ir constant

M=[];
for k=1:length(ParNames)
    
    LineIndex=inpFindString(data,ParNames{k},'printerror','no');
    
    if isnan(LineIndex); M(k,1)=NaN; continue; end
    
    charIndex=strfind(data{LineIndex},':');
    M(k,1)=inpReadNumber(data,LineIndex,'%f','skip',charIndex);

    M(k,1)=M(k,1)/ScaleFactor(k);
    
end

% % M_out = struct('names',parNamesRed);
% [M_out(:).data] = deal(M);

M_out=struct();

for k=1:length(ParNames)
	M_out=setfield(M_out,ParNamesRed{k},M(k));
end





