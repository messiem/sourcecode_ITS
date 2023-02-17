function [output,age_structure]=its_agestructured_model(forcing_timeseries,m)


%% ITS_AGESTRUCTURED_MODEL: computes the age structure of a population characterized by a natural mortality rate of m,
%% and force it using forcing_timeseries as an input where animals are born proportional to forcing_timeseries when positive
% The age structure is defined with daily resolution, and forcing_timeseries is similarly assumed to be daily.
%
% [output,age_structure]=its_agestructured_model(forcing_timeseries,m)
%
% INPUTS:
% forcing_timeseries 	daily forcing (only positive values are considered)
% m						natural mortality rate, in year^-1
%
% OUTPUTS:
% output				structure containing
%							.counts: numbers for each age group and time step (age x time, both with daily resolution)
%							.counts_all: sum over all age groups
% age_structure			structure defining the model as parameterized by m, containing
%							.death_rate: probability of dying at each age
%							.age: age of each age group, in days, starting at 0
%							.age_pyramid: survivorship function, ie. proportion alive at each age (starting at 1 for age 0)
%							.mean_LE: mean life expectancy at birth
%
% Monique Messié, September 2022 for public version


% compute age structure
age_structure=struct();
age_structure.death_rate=1-exp(-m/365.25);				% proportion dying every day
age_structure.age=(0:round(100*365.25))';
age_structure.age_pyramid=zeros(length(age_structure.age),1);
age_structure.age_pyramid(1)=1;
for iage=2:length(age_structure.age)
	age_structure.age_pyramid(iage)=(1-age_structure.death_rate).*age_structure.age_pyramid(iage-1);	% all individuals age, except for those who die
	if age_structure.age_pyramid(iage)<=0.01			% lifespan: when numbers drop below 1% of natality
		age_structure.lifespan=age_structure.age(iage); 
		age_structure.age=age_structure.age(1:iage);	% stop calculations after lifespan (i.e., assume everyone dies then)
		age_structure.age_pyramid=age_structure.age_pyramid(1:iage);	
		break
	end
end


% Compute mean life expectancy following https://www.lifeexpectancy.org/lifetable.shtml
death=[age_structure.age_pyramid(1:end-1)-age_structure.age_pyramid(2:end);age_structure.age_pyramid(end)];	% number of deaths from age x to x+1
L=[age_structure.age_pyramid(2:end)+0.5*death(1:end-1);0];	% total number of person-years lived by the cohort from age x to x+1
age_structure.mean_LE=sum(L);								% mean life expectancy at birth


% Force the model using the forcing_timeseries (only when positive)
forcing_timeseries(forcing_timeseries<0)=0;
output=struct();
output.counts=zeros(length(age_structure.age),length(forcing_timeseries));
output.counts(:,1)=age_structure.age_pyramid*mean(forcing_timeseries);							% starting with mean age pyramid
for itime=2:length(forcing_timeseries)
	output.counts(1,itime)=forcing_timeseries(itime);											% natality @age 1
	output.counts(2:end,itime)=(1-age_structure.death_rate).*output.counts(1:end-1,itime-1);	% all individuals age, except for those who die
end
output.total_counts=squeeze(sum(output.counts,1))';


return


