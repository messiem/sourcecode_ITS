function varargout=its_compute_PCA(taxa,varargin)


%% ITS_COMPUTE_PCA: for a given taxa structure as obtained from its_load_data, calculate PCA, 
%%  returning the first PC and corresponding loadings.
%
% [PCA_mode1,PC1_monthly]=its_compute_PCA(taxa,varargin)
%
% REQUIRED INPUT:
% taxa		 	structure obtained from its_load_data, with additional field .counts_norm 
%				(typically 4th-root and mean removed, see start_sourcecode_ITS)
%
% OPTIONAL INPUT:
% 'plot'		displays the time series as in Fig. 1
%
% OUTPUT:
% PCA_mode1		structure with PC1 mode 1 PCA_mode1s: fields time, PC1, mode1_loadings, and varex1
% PC1_monthly	structure containing PCA_mode1.PC1 regridded monthly with 1-gap filled: fields time, PC1
%
% Monique Messié, September 2022 for public version


if ~isfield(taxa,'counts_norm')
	error('Provide counts_norm within the input structure, such as 4th-root transformed and mean removed')
end


% Compute PCA on normalized counts (removing the long-term mean)
[u,s,v]=svd(taxa.counts_norm-repmat(mean(taxa.counts_norm,2),1,length(taxa.time)),'econ'); 
egn=diag(s*s');       						% eigenvalues
loadings=u';								% eigenfunctions
PCs=v.*repmat(sqrt(egn)',size(v,1),1); 		% principal component
varex=egn./sum(egn)*100; 					% variance explained


% Only keep mode 1
PCA_mode1=struct();
PCA_mode1.time=taxa.time;
PCA_mode1.PC1=PCs(:,1);
PCA_mode1.mode1_loadings=loadings(1,:);
PCA_mode1.varex1=varex(1);
if mean(PCA_mode1.mode1_loadings)<0			% since the sign is arbitrary, flip it if necessary to have positive loadings
	PCA_mode1.mode1_loadings=-PCA_mode1.mode1_loadings;
	PCA_mode1.PC1=-PCA_mode1.PC1;
end
norm_factor=max(abs(PCA_mode1.PC1)); 			% normalizing PC1 (and correcting mode1_loadings as well)
PCA_mode1.PC1=PCA_mode1.PC1/norm_factor; 
PCA_mode1.mode1_loadings=PCA_mode1.mode1_loadings*norm_factor;


% Compute monthly PCA_mode1, used for calulating p-values adjusted for autocorrelation and in figures
PC1_monthly=struct();
[y,m,~]=datevec(PCA_mode1.time);
PC1_monthly.time=(datetime(y(1),m(1):m(1)+12*(y(end)-y(1)+1),15))';
PC1_monthly.time=PC1_monthly.time(PC1_monthly.time<=datetime(y(end),m(end),15));
PC1_monthly.PC1=nan(size(PC1_monthly.time));
for itime_monthly=1:length(PC1_monthly.time)
	itime_daily=datetime(y,m,15)==PC1_monthly.time(itime_monthly);		% identify daily timesteps corresponding to each month
	PC1_monthly.PC1(itime_monthly)=mean(PCA_mode1.PC1(itime_daily));			% average corresponding data
end


% Fill 1-month gaps in PC1_monthly
list_gaps=find(isnan(PC1_monthly.PC1));
list_gaps=list_gaps(list_gaps>1 & list_gaps<length(PC1_monthly.PC1));
for igap=list_gaps'
	PC1_monthly.PC1(igap)=mean([PC1_monthly.PC1(igap-1);PC1_monthly.PC1(igap+1)]);	%will remain at NaN if value before or after is NaN
end


% Display the result (Fig. 2)
if max(strcmp(varargin,'plot'))

	% Sort loadings
	[mode1_loadings,iok]=sort(PCA_mode1.mode1_loadings); 
	taxaname=taxa.taxaname(iok); 
	groupname=taxa.groupname(iok); 

	% Set up figure
	group_list=unique(taxa.groupname);
	color_list=jet(length(group_list));
	h=nan(length(group_list),1);				% keep track of handles to use for the legend
	figure
	gf=gcf;
	gf.Position(3)=gf.Position(3)*0.7;

	% Time series (monthly 1-gap-filled PC1)
	axes('Position',[0.1 0.83 0.8 0.12],'FontSize',8), hold on
		if strcmp(taxa.community,'benthos')
			plot(PC1_monthly.time,PC1_monthly.PC1,'Marker','*','LineStyle','none')
		else
			plot(PC1_monthly.time,PC1_monthly.PC1)
		end
		plot(xlim,[0 0],'k')
		ylim([-1 1])
		title([taxa.community,' mode 1 (',num2str(round(PCA_mode1.varex1)),'%)'],'FontSize',10)

	% Loadings
	axes('Position',[0.45 0.04 0.4 0.73]), hold on
		for itaxa=1:length(taxaname)
			igroup=strcmp(group_list,groupname{itaxa});
			rgbcolor=color_list(igroup,:); 
			h(igroup)=patch([0,0,mode1_loadings(itaxa),mode1_loadings(itaxa)],...
				[itaxa-0.4,itaxa+0.4,itaxa+0.4,itaxa-0.4],rgbcolor);
			text(min(mode1_loadings)*1.1,itaxa,taxaname{itaxa},'FontSize',6,'HorizontalAlignment','right','Color',rgbcolor)
		end
		set(gca,'YLim',[0.5 length(taxaname)+0.5],'YTick',1:length(taxaname),'YTickLabel','','YDir','reverse','FontSize',8), grid on
	legend(h,group_list,'FontSize',6)

end


% Return the result
varargout={PCA_mode1,PC1_monthly};
varargout=varargout(1:nargout);

return
