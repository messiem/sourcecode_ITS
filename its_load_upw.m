function upw=its_load_upw(windstation)


%% ITS_LOAD_UPW: loads the upwelling time series for a wind station
%
% upw=its_load_upw(windstation);
%
% REQUIRED INPUT:
% windstation 	'MontereyBay' (for surface and midwater) or 'PtConception' (for benthos)
%
% OUTPUT:
% upw		structure with fields windstation, time (Matlab format), upw (upwelling), unit
%
% Monique Messié, September 2022 for public version



% Load upwelling data
filename=['data/upwelling_',windstation,'.csv'];
data_csv=readcell(filename);
headers_all=data_csv(1,:);
data_all=cell2mat(data_csv(2:end,:));


% Put data into a structure
upw=struct();
upw.windstation=windstation;
upw.time=datetime(data_all(:,1),data_all(:,2),data_all(:,3));	% headers_all(2,1:3) are year, month, day
upw.upw=data_all(:,4);				% formatted as upw x time
upw_header=headers_all{4};
upw.unit=upw_header(strfind(upw_header,'(')+1:strfind(upw_header,')')-1);


return