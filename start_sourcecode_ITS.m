%% START_SOURCECODE_ITS: code to reproduce figures and results from the paper
% These functions must be run within the sourcecode_ITS folder for data to load correctly.
%
% Reference: Messié, M., R.E. Sherlock, C.L. Huffard, J.T. Pennington, C.A. Choy, R.P. Michisaki, K. Gomes, F.P. Chavez, B.H. Robison, and K.L. Smith Jr. (2023)
% 	Coastal upwelling drives ecosystem temporal variability from the surface to the abyssal seafloor, Proceedings of the National Academy of Sciences, in press.


community='midwater';		% set to surface, midwater, or benthos


% Determine if Statistics and Machine Learning Toolbox is available (required to compute p-values)
V = ver;
VName = {V.Name};
if any(strcmp(VName, 'Statistics and Machine Learning Toolbox')), compute_pvalue=true; 
else, compute_pvalue=false; 
	disp('p-values won''t be computed; requires the Statistics Toolbox')
end


% Load data and reproduce Fig. 1
taxa=its_load_taxa(community,'plot');
print('-djpeg','-r300',['outputs/Fig1_taxa_',community,'.jpg'])


% Data pre-processing as explained in Messié et al.
if strcmp(taxa.community,'benthos')			% restrict benthic time series to 1996-present to ensure PCA mode separation
	iok=taxa.time>=datetime(1996,1,1);
	taxa.time=taxa.time(iok);
	taxa.counts=taxa.counts(:,iok);
end	
taxa.counts_norm=taxa.counts.^(1/4);		% 4th-root transform


% Compute PCA and reproduce Fig. 2 (monthly PCA is used for calulating p-values adjusted for autocorrelation and in figures)
[PCA_mode1,PC1_monthly]=its_compute_PCA(taxa,'plot');
print('-djpeg','-r300',['outputs/Fig2_PCA_mode1_',community,'.jpg'])


% Load daily upwelling time series, compute monthly time series, and identify time steps coincident with the taxa time series
% For the benthos time series that is on a monthly time scale, the coincident time steps are based on the monthly upwelling time series.
switch community
case {'surface','midwater'}
	upw=its_load_upw('MontereyBay');					% upwelling time series
	[y,m,~]=datevec(upw.time);
	time_monthly=unique(datetime(y,m,15));				
	itime_taxa=ismember(upw.time,PCA_mode1.time);		% upwelling time steps coinciding with taxa time series data
case 'benthos'
	upw=its_load_upw('PtConception');
	[y,m,~]=datevec(upw.time);
	time_monthly=unique(datetime(y,m,15));
	itime_taxa=ismember(time_monthly,PCA_mode1.time);	% (taxa time is monthly, set on the 15th of the month)
end


% Compute upwelling integration curve
% Setting tau resolution to daily until 3 weeks, then weekly until 3 months, then monthly until 4.5 years, then yearly.
corr_intupw=struct('tau',[1:21 28:7:90 120:30.5:365.25*4.5 (4.5:0.5:10)*365.25]'); % in days (upwelling is a daily time series)
corr_intupw.r=nan(size(corr_intupw.tau));
corr_intupw.p=nan(size(corr_intupw.tau));
upw_anom=upw.upw-mean(upw.upw);
for itau=1:length(corr_intupw.tau)
	intupw=its_compute_integration(upw_anom,corr_intupw.tau(itau));
	intupw_monthly=nan(size(time_monthly));
	for itime_monthly=1:length(time_monthly)
		itime_daily=datetime(y,m,15)==time_monthly(itime_monthly);	% identify daily timesteps corresponding to each month
		intupw_monthly(itime_monthly)=mean(intupw(itime_daily));	% average corresponding data
	end
	if strcmp(community,'benthos'), intupw=intupw_monthly; end		% original time resolution is monthly for the benthos
	% Compute correlation using original time resolution
	r=corrcoef(intupw(itime_taxa),PCA_mode1.PC1);
	corr_intupw.r(itau)=r(1,2); 
	% Compute p-value ajusted for autocorrelation using monthly, 1-gap filled time series (regular time steps required)
	if compute_pvalue
		[~,p]=its_compute_adjusted_correlation(intupw_monthly(ismember(time_monthly,PC1_monthly.time)),...
												PC1_monthly.PC1(ismember(PC1_monthly.time,time_monthly)));
		corr_intupw.p(itau)=p; 
	end
end


% Identify damping timescale and compute corresponding integration 
% (for midwater, use tau = 4.5 years since the damping timescale is undefined)
tau_damping = min(corr_intupw.tau(corr_intupw.r==max(corr_intupw.r)));
if tau_damping==max(corr_intupw.tau), tau_damping=NaN; end			% asymptomatic curve so undefined
disp(['Damping timescale: ',num2str(tau_damping),' days (',num2str(tau_damping/365.25),' years)'])
if isnan(tau_damping) && strcmp(community,'midwater'), tau=4.5*365.25; else, tau=tau_damping; end
intupw=its_compute_integration(upw_anom,tau);


% Reproduce Fig. 3
figure
subplot(2,1,1), hold on
	h1=plot(corr_intupw.tau/365.25,corr_intupw.r,'b:','LineWidth',1.5);
	h2=plot(corr_intupw.tau(corr_intupw.p<0.01)/365.25,corr_intupw.r(corr_intupw.p<0.01),'b','LineWidth',1.5);
	ylim([0 1])
	if strcmp(community,'surface'), xlim([0 2]), end		% zooming in to better see surface results
	h3=plot([1 1]*tau_damping/365.25,ylim,'b--','LineWidth',2);
	xlabel('Upwelling integration timescale (years)','FontSize',10)
	ylabel('Correlation PC1 vs intupw','FontSize',10)
	legend([h1,h2,h3],{'correlation (p<0.01)',['correlation (p' char(8805) '0.01)'],'\tau_{damping}'},...
		'Location',[0.7 0.65 0.1 0.05])
	title(['Integration curve ',community])
subplot(2,1,2), hold on
	switch community
	case 'benthos', h1=plot(PCA_mode1.time,PCA_mode1.PC1,'Marker','*','LineStyle','none');
	otherwise, h1=plot(PCA_mode1.time,PCA_mode1.PC1,'Marker','*');
	end
	h2=plot(upw.time,intupw/max(abs(intupw)),'LineWidth',1);
	xlim([min(upw.time) max(upw.time)])
	plot(xlim,[0 0],'k')
	ylim([-1 1])
	xlabel('Time','FontSize',10)
	ylabel('normalized','FontSize',10)
	legend([h1,h2],{[community,' PC1'],...
					['Integrated upwelling (r=',num2str(round(corr_intupw.r(corr_intupw.tau==tau),2)),')']},...
			'Location',[0.3 0.38 0.1 0.05])
	title('PC1 and integrated upwelling time series')
print('-djpeg','-r300',['outputs/Fig3_PC1_vs_intupw_',community,'.jpg'])


% Find best model parameterization that reproduces a taxa time series (reproduce Fig. 5 top & middle panels)
if strcmp(community,'midwater')

	for taxaname={'Aegina+','Cyclothone'}, taxaname=taxaname{:};

		% Define best model parameterization
		% (for other taxa, the best model can be found by looping through various values of m and finding the one that maximizes r_model_taxa below)
		switch taxaname
		case 'Aegina+'; m=1.3;
		case 'Cyclothone', m=0.25;
		otherwise, error('Best model parameterization is unknown, find it and update the code')
		end

		% Extract the taxonomic time series
		taxa_series=taxa.counts_norm(strcmp(taxa.taxaname,taxaname),:);

		% Compute modeled time series and correlation with taxa_series 
		[output,age_structure]=its_agestructured_model(upw.upw,m);
		output.time=upw.time;
		taxa_model=output.total_counts.^(1/4);
		r=corrcoef(taxa_model(itime_taxa),taxa_series); 
		r_model_taxa=r(1,2);

		% Compute the integration curves for measured and modeled taxonomic time series, following the same method as above
		corr_intupw=struct('tau',[1:21 28:7:90 120:30.5:365.25*4.5 (4.5:0.5:10)*365.25]'); % in days (upwelling is a daily time series)
		corr_intupw.r_taxa=nan(size(corr_intupw.tau));
		corr_intupw.r_model=nan(size(corr_intupw.tau));
		upw_anom=upw.upw-mean(upw.upw);
		iperiod=upw.time>=taxa.time(1) & upw.time<=taxa.time(end);	% compute the integration curve based on the same time period as data
		for itau=1:length(corr_intupw.tau)
			intupw=its_compute_integration(upw_anom,corr_intupw.tau(itau));
			r=corrcoef(intupw(itime_taxa),taxa_series);
			corr_intupw.r_taxa(itau)=r(1,2); 
			r=corrcoef(intupw(iperiod),taxa_model(iperiod));
			corr_intupw.r_model(itau)=r(1,2); 
		end

		% Find taxa-specific damping timescale (and same for the model) and compute the corresponding wind integration
		tau_damping_taxa = min(corr_intupw.tau(corr_intupw.r_taxa==max(corr_intupw.r_taxa)));
		tau_damping_model = min(corr_intupw.tau(corr_intupw.r_model==max(corr_intupw.r_model)));
		intupw=its_compute_integration(upw_anom,tau_damping_taxa);

		% Figure set up
		ydata=mean(taxa_series)+[-1 1]*max(abs(taxa_series-mean(taxa_series)));		% compute ylims such that values are centered on the series' mean
		ymodel=mean(taxa_model)+[-1 1]*max(abs(taxa_model-mean(taxa_model)));
		figure

		% Top panel: time series
		xlims=[min(upw.time) max(upw.time)];
		ax1=axes('Position',[0.1 0.6 0.8 0.37],'YAxisLocation','left');
		h1=line(taxa.time,taxa_series,'Color','r','LineWidth',2);
		set(gca,'XLim',xlims,'YLim',ydata)
		ylabel(taxaname)
		axes('Position',ax1.Position,'Color','none','YAxisLocation','right')
		h2=line(upw.time,intupw,'Color','k','LineWidth',0.5);
		set(gca,'XTick',[],'XLim',xlims,'YLim',[-2 2])
		ylabel('Integrated upwelling')
		axes('Position',ax1.Position,'Color','none')
		h3=line(output.time,taxa_model,'Color',[0 0.7 1],'LineWidth',0.5);
		set(gca,'XTick',[],'XLim',xlims,'YLim',ymodel,'YTick',[])
		legend([h1,h2,h3],{[taxaname,' time series'],...
							['Integrated upwelling (r=',num2str(round(max(corr_intupw.r_taxa),2)),')'],...
							['Modeled time series (r=',num2str(round(r_model_taxa,2)),')']},...
				'Location',[0.25 0.92 0.05 0.05])

		% Bottom panel: integration curves
		axes('Position',[0.1 0.1 0.8 0.37]), hold on
		h1=plot(corr_intupw.tau/365.25,corr_intupw.r_taxa,'Color','k','LineWidth',2);
		h2=plot(corr_intupw.tau/365.25,corr_intupw.r_model,'Color',[0 0.7 1],'LineWidth',2);
		ylim([0 1])
		plot([1 1]*tau_damping_taxa/365.25,ylim,'k--','LineWidth',2)
		plot([1 1]*tau_damping_model/365.25,ylim,'--','Color',[0 0.7 1],'LineWidth',1)
		xlabel('Upwelling integration timescale (years)')
		ylabel('Integration curve')
		legend([h1,h2],{['data (',taxaname,')'],'model'},'Location',[0.75 0.15 0.05 0.05])
		print('-djpeg','-r300',['outputs/Fig5_model_vs_intupw_vs_',taxaname,'.jpg'])

	end

end