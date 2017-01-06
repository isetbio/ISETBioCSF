function [fractionCorrect,dPrime] = analyticPoissonIdealObserver(alphaMeanResponses,betaMeanResponses)
%analyticPoissonIdealObserver]  Get the ideal observer fraction correct in a TAFC task
%    fractionCorrect = analyticPoissonIdealObserver(alphaMeanResponses,betaMeanResponses)
%
%    Get the fraction correct in a TAFC task using the formula for a
%    Poisson ideal observer developed by Geisler, 1984, JOSA A, 1, pp. 775
%    ff.
%
%    Also returns dPrime.  Indeed, the Geisler formula returns dPrime and
%    this routine then converts to TAFC percent correct by numerically
%    integrating the area under the ROC curve for that dPrime, assuming
%    normal distributions
%     with equal variance in the conversion.
%
%    The two input arguments should all be vectors of the same length, and
%    give the mean of the Poisson distributed responses for the two
%    stimulus types.
%
% See also dPrimeToTAFCFractionCorrect

% Handle special case when the mean responses for the two classes are
% identical, which runs into numerical problems if we compute it using the
% code below
if (all(alphaMeanResponses == betaMeanResponses))
    fractionCorrect = 0.5;
    dPrime = 0;
    return;
end

% This comes from the appendix of the paper
numerator = sum( (betaMeanResponses-alphaMeanResponses).*log(betaMeanResponses./alphaMeanResponses) );
denominator = 0.5*sum( (betaMeanResponses+alphaMeanResponses).*(log(betaMeanResponses./alphaMeanResponses).^2) );
dPrime = numerator / sqrt(denominator);

% Call into our dPrime to TAFC conversion function
fractionCorrect = dPrimeToTAFCFractionCorrect(dPrime);



