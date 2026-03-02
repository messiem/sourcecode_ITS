function integrated_series=its_compute_integration(input_series,tau)


%% ITS_COMPUTE_INTEGRATION: integrates a time series with a given integration timescale following Di Lorenzo and Ohman (2013).
% Reference: Di Lorenzo, E., & Ohman, M.D. (2013). A double-integration hypothesis to explain ocean ecosystem response to climate forcing. 
% Proceedings of the National Academy of Sciences, 110(7), 2496-2499, https://doi.org/10.1073/pnas.1218022110
%
% IMPORTANT: the integration timescale tau is in units of time steps.
% So if the input timeseries is daily, tau should be expressed in days; 
% if monthly, it should be expressed in months, etc.
%
% integrated_series=its_compute_integration(input_series,tau)
%
% REQUIRED INPUTS:
% input_series		input time series on regular time steps, with no gaps
% tau				integration time scale, expressed in time steps
%
% OUTPUT:
% integrated_series		input_series integrated with timescale tau
%
% Monique Messié, September 2022 for public version
% January 2026: modified equation so that integrated_series(ipts) is a function of input_series(ipts) rather than input_series(ipts-1), as in other AR(1) equations
%					(see eg https://doi.org/10.1038/s42003-021-02960-y)
%				This doesn't change anything to the paper results where integration time scales were of several days at least.
%				Also modified standardization so that the first data point remains unchanged (no impact to results based on correlations).


% check that the input time series is gap-free
if max(isnan(input_series)), error('need a gap-free input_series for integration!!'), end


% compute integration
integrated_series=input_series*NaN; 
integrated_series(1)=input_series(1);
for ipts=2:length(input_series)
	integrated_series(ipts)=input_series(ipts)+integrated_series(ipts-1)*(1-1/tau);
end


% standardize the standard deviation to the original input_series standard deviation, ensuring the first data point remains identical (cf l. 32)
integrated_series=input_series(1)+(integrated_series-input_series(1))/std(integrated_series)*std(input_series);


return