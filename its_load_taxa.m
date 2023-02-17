function taxa=its_load_taxa(community,varargin)


%% ITS_LOAD_TAXA: for a given community (surface, midwater, or benthos), returns the time series in a structure.
%
% taxa=its_load_taxa(community,varargin);
%
% REQUIRED INPUT:
% community 	'surface','midwater', or 'benthos'
%
% OPTIONAL INPUT:
% 'plot'		displays the time series as in Fig. 1
%
% OUTPUT:
% taxa		structure with fields community, time (Matlab format), taxaname (taxonomic group name), 
%				groupname (higher-level group name), counts (data formatted as taxa x time), unit
%
% Monique Messié, September 2022 for public version



% Load time series data
filename=['data/',community,'.csv'];
data_csv=readcell(filename);
headers_all=data_csv(1:2,:);
data_all=cell2mat(data_csv(3:end,:));


% Put data into a structure
taxa=struct();
taxa.community=community;
switch community
case 'benthos'	% monthly data: the file only provides year and month (setting day to the middle of the month)
	taxa.time=datetime(data_all(:,1),data_all(:,2),15);				% headers_all(2,1:2) are year, month
	first_col=3;													% first column that contains a taxon time series
otherwise
	taxa.time=datetime(data_all(:,1),data_all(:,2),data_all(:,3));	% headers_all(2,1:3) are year, month, day
	first_col=4;
end
taxa.taxaname=headers_all(2,first_col:end);
taxa.groupname=headers_all(1,first_col:end);
taxa.counts=data_all(:,first_col:end)';				% formatted as taxa x time
taxa.unit=headers_all(1);


% Go back to original unit if needed (midwater density)
if contains(taxa.unit,'x1000')
	taxa.counts=taxa.counts/1000; 
	taxa.unit=strrep(taxa.unit,' (x1000)','');
end


% Plot as in Fig. 1
% The colorbar in Fig. 1 is based on 'balance' from cmocean (https://doi.org/10.5670/oceanog.2016.66) instead of 'turbo'.
if max(strcmp(varargin,'plot'))

	counts_norm=taxa.counts.^(1/4);												% 4th-root transform
	counts_norm=counts_norm-repmat(mean(counts_norm,2),1,length(taxa.time));	% remove long-term mean
	[time2D,name2D]=meshgrid(taxa.time,1:length(taxa.taxaname));
	switch community
	case 'surface', dataLims=[-3 3];
	case 'midwater', dataLims=[-0.25 0.25];
	case 'benthos', dataLims=[-0.8 0.8];
	end

	figure
	scatter(time2D(:),name2D(:),6,counts_norm(:),'filled')
	caxis(dataLims)	
	shading flat
	colorbar
	colormap('turbo')
	ylim([0 length(taxa.taxaname)+1])
	set(gca,'YTick',1:length(taxa.taxaname),'YTickLabel',taxa.taxaname,'YDir','reverse','FontSize',6)
	xlabel('Time','FontSize',8)
	title([community,' time series'],'FontSize',10)
	
end


return