% A criterion：MFE saturation level = function of temperature  
% if we assume saturation criteria is linear function of temperature, then we can write the fuction as,

function MFEsat = mfesat(T,cal)
% MFESAT  Compute MFE saturation level as a function of temperature
%
% Usage:
%   MFEsat = mfesat(T,cal)
%
% Input:
%       T            - temperature (°C, scalar or vector)
%   cal.MFE_liquidus - saturation value at liquidus (wt fraction)
%   cal.T_liquidus   - liquidus temperature (°C)
%   cal.T_solidus    - solidus temperature (°C)
%   cal.k            - empirical coefficient (0–1)
%
% Output:
%   MFEsat           - saturation level (wt fraction)
%

% ---- Calculate temperature-dependent saturation ----
MFEsat = cal.MFE_liquidus .* ...
    (1 - cal.k .* (cal.T_liquidus - T) ./ ...
    (cal.T_liquidus - cal.T_solidus));

% ---- Constrain within valid range ----
MFEsat = max(MFEsat, 1e-4);                  % no negative values
MFEsat = min(MFEsat, cal.MFE_liquidus);      % cannot exceed liquidus value
end

