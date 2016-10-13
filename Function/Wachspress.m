function [phi, dphi] = Wachspress(v, x)
% [phi, dphi] = Wachspress(v, x)
%==========================================================================
% Reference: GRADIENT BOUNDS FOR WACHSPRESS COORDINATES ON POLYTOPES
%            MICHAEL S. FLOATER, ANDREW GILLETTE, N. SUKUMAR
%
% Evaluate Wachspress basis functions in a convex polygon
%
% Inputs:
% v : [x1 y1; x2 y2; ...; xn yn], the n vertices of the polygon in ccw
% x : [x(1) x(2)], the point at which the basis functions are computed
% Outputs:
% phi : output basis functions = [phi_1; ...; phi_n]
%==========================================================================

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

n = size(v, 1);
w = zeros(n, 1);
%phi = zeros(n,1);
dphi = zeros(n, 2);
p = zeros(n, 2);
for i = 1 : n
    d = v(mod(i, n) + 1, :) - v(i, :);
    normal_vector = [d(2), -d(1)] / norm(d);
    h = (v(i, :) - x) * normal_vector'; % dot product
    p(i, :) = normal_vector / h;  % scaled normal vectors
end
R = zeros(n, 2);
for i = 1 : n
    im1 = mod(i - 2, n) + 1;
    w(i) = det([p(im1, :); p(i, :)]);
    R(i, :) = p(im1, :) + p(i, :);
end
wsum = sum(w);
phi = w / wsum;

phiR = phi' * R;
for k = 1 : 2
    dphi(:, k) = phi .* (R(:, k) - phiR(:, k));
end
end
