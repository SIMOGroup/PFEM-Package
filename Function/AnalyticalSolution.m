function [eldisp, estress, estrain] = AnalyticalSolution(x, y, p, emodu, poisson, d, l)
% Example: A Rectangular Cantilever Loaded at the end
%  ux=P*y/(6EI)*[(6L-3x)*x+(2+poisson)*(y^2-D^2/4)];
%  uy=-P/(6EI)*[3*poisson*y^2*(L-x)+(4+5*poisson)*D^2*x/4+(3*L-x)*x^2];
%  sigmaX=P*(L-x)*y/I
%  sigmaY=0
%  sigmaXY=P*(y*y-D*D/4)/(2*I);
%--------------------------------------------------------------------------
%{
Copyright (C) <2016>  <Son Nguyen-Hoang, Khai Chau-Nguyen, Hung Nguyen-Xuan>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

I = d*d*d / 12;
PEI = p / (6*emodu*I);
ux = PEI*y*((6*l - 3*x)*x + (2 + poisson)*(y*y - d*d / 4));
uy = -PEI*(3*poisson*y*y*(l - x) + (4 + 5*poisson)*d*d*x / 4 + (3*l - x)*x*x);
dux = PEI*6*(l - x)*y;
dvy = -PEI*6*poisson*y*(l - x);
dxy = PEI*6*(1 + poisson)*(y*y - d*d / 4);               
sigx = p*(l - x)*y / I;
sigy = 0;
sigxy = p*(y*y - d*d / 4) / 2 / I;
eldisp = [ux; uy];
estrain = [dux; dvy; dxy];
estress = [sigx; sigy; sigxy];
end

