function varargout=its_compute_adjusted_correlation(x1,x2)


%% ITS_COMPUTE_ADJUSTED_CORRELATION: computes Pearson correlations between 2 time series 
%% and the p-value adjusted for autocorrelation following Pyper & Peterman (1998).
% Reference: Pyper, B.J., & Peterman, R.M. (1998). Comparison of methods to account for autocorrelation in correlation analyses of fish data. 
% 			Canadian Journal of Fisheries and Aquatic Sciences, 55(9), 2127-2140, https://doi.org/10.1139/f98-104 
% The method follows the paper recommendation: "Specifically, we recommend that researchers use eq. 1 without the weighting function, 
%			with autocorrelations estimated over N/5 lags j using eq. 7, and with critical value using N* â€“ 2 degrees of freedom."
%
% [r,adjusted_p]=its_compute_adjusted_correlation(x1,x2);
%
% REQUIRED INPUTS:
% 	x1, x2 vectors of the same length ON REGULAR TIME STEPS (e.g., monthly here).
%
% OUTPUTS:
% 	r	Pearson correlation coefficient
% 	p	p-value adjusted for autocorrelation following Pyper & Peterman 1998
%
% Monique Messié, December 2022


% Ensure time series are on the same time period (different gaps are OK, but autocorrelation should be at least estimated over the same time period)
idata = ~isnan(x1) | ~isnan(x2);
iperiod = find(idata,1,'first'):find(idata,1,'last');
x1=x1(iperiod); x2=x2(iperiod);


% Normalise time series to simplify autocorrelation equations
x1norm = (x1-mean(x1,'omitnan'))/std(x1,'omitnan');
x2norm = (x2-mean(x2,'omitnan'))/std(x2,'omitnan');


% Compute autocorrelation using eqn 7 in Pyper & Peterman (modified Chelton)
N = length(x1); 
nlags=floor(N/5);
autocorr=struct();
autocorr.x1=nan(nlags,1);
autocorr.x2=nan(nlags,1);
for j=1:nlags
	autocorr.x1(j) = N/(N-j)/(N-1)*sum(x1norm(1:N-j).*x1norm(j+1:N),'omitnan');		
	autocorr.x2(j) = N/(N-j)/(N-1)*sum(x2norm(1:N-j).*x2norm(j+1:N),'omitnan');	
end
autocorr.x1(autocorr.x1==0)=NaN;		% because sum([NaN NaN],'omitnan')=0
autocorr.x2(autocorr.x2==0)=NaN;


% Compute Pearson correlation coefficient
iok = ~isnan(x1) & ~isnan(x2);
if sum(iok)>1
	r=corrcoef(x1(iok),x2(iok));
	r=r(1,2);
else, r=NaN;
end


% Compute DF (degrees of freedom) using eqn 1 in Pyper & Peterman 
Ndata = sum(~isnan(x1) & ~isnan(x2));
DF = Ndata/(1+2*sum(autocorr.x1.*autocorr.x2,'omitnan'));	
if DF>Ndata, DF=Ndata; end


% Computed adjusted p-value using adjusted degrees of freedom = DF-2
adjusted_p=(1-fcdf((1+abs(r))/(1-abs(r)),DF-2,DF-2))*2;	


% Return the result
varargout={r,adjusted_p};
varargout=varargout(1:nargout);


end

