function [weight, point] = TriangularQuad(p, degreeGauss)
% [weight, point] = TriangularQuad(p, degreeGauss)
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

[Q, W] = Quadrature(1, degreeGauss);
point = zeros(length(W), 2); 
weight = zeros(length(W), 1);
for q = 1 : length(W)
    xi = Q(q, 1); 
    et = Q(q, 2);
    N = [1 - (xi + et), xi, et];
    Nxi = [-ones(size(xi)), ones(size(xi)), zeros(size(xi))];
    Net = [-ones(size(xi)), zeros(size(xi)), ones(size(xi)),];
    dNds = [Nxi; Net];
    J0 = p'*dNds';
    point(q, :) = N*p;
    weight(q) = det(J0)*W(q);
end
end