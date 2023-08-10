function  res = TV_2DPERF()

%res = TVOP()
%
% Implements a spatial finite-differencing operator for dynamic data.
%
% Ricardo Otazo 2008

res.adjoint = 0;
res = class(res,'TV_2DPERF');

