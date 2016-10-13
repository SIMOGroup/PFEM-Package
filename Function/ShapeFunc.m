function [N, Nxi, Neta] = ShapeFunc(elem, nen, pospg) 
% [N, Nxi, Neta] = ShapeFunc(elem, nen, pospg)
%
% Input:    
%   elem:   Type of element (0 for quadrilaterals and 1 for triangles)
%   nen:    number of element nodes
%   pospg:  coordinates of Gauss points in the reference element
%
% Output:
%   N, Nxi, Neta: matrices storing the values of the shape functions on the Gauss points
%                 of the reference element. Each row concerns to a Gauss point
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

xi = pospg(:, 1); 
et = pospg(:, 2); 

if elem == 0   % quadrilaterals
    if nen == 1  %P0
        N    = ones(size(xi));
        Nxi  = zeros(size(xi));
        Neta = zeros(size(xi));
    elseif nen == 4  %Q1
        N    = [(1-xi).*(1-et)/4, (1+xi).*(1-et)/4, ...
               (1+xi).*(1+et)/4, (1-xi).*(1+et)/4]; 
        Nxi  = [(et-1)/4, (1-et)/4, (1+et)/4, -(1+et)/4]; 
        Neta = [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ]; 
    elseif nen == 9  %Q2
        N    = [xi.*(xi-1).*et.*(et-1)/4, xi.*(xi+1).*et.*(et-1)/4, ...
               xi.*(xi+1).*et.*(et+1)/4, xi.*(xi-1).*et.*(et+1)/4, ...
               (1-xi.^2).*et.*(et-1)/2,  xi.*(xi+1).*(1-et.^2)/2,   ...
               (1-xi.^2).*et.*(et+1)/2,  xi.*(xi-1).*(1-et.^2)/2,   ...
               (1-xi.^2).*(1-et.^2)];
        Nxi  = [(xi-1/2).*et.*(et-1)/2,   (xi+1/2).*et.*(et-1)/2, ...
               (xi+1/2).*et.*(et+1)/2,   (xi-1/2).*et.*(et+1)/2, ...
               -xi.*et.*(et-1),          (xi+1/2).*(1-et.^2),   ...
               -xi.*et.*(et+1),          (xi-1/2).*(1-et.^2),   ...
               -2*xi.*(1-et.^2)];
        Neta = [xi.*(xi-1).*(et-1/2)/2,    xi.*(xi+1).*(et-1/2)/2, ...
               xi.*(xi+1).*(et+1/2)/2,    xi.*(xi-1).*(et+1/2)/2, ...
               (1-xi.^2).*(et-1/2),       xi.*(xi+1).*(-et),   ...
               (1-xi.^2).*(et+1/2),       xi.*(xi-1).*(-et),   ...
               (1-xi.^2).*(-2*et)];
   else
       error ('Error in ShapeFunc: unavailable quadrilateral')
   end
elseif elem == 1  % triangles
    if nen == 1  %P0
        N    = ones(size(xi));
        Nxi  = zeros(size(xi));
        Neta = zeros(size(xi));
    elseif nen == 3  %P1
        N    = [1-(xi+et),xi,et]; 
        Nxi  = [-ones(size(xi)),ones(size(xi)),zeros(size(xi))]; 
        Neta = [-ones(size(xi)),zeros(size(xi)),ones(size(xi)),]; 
        
    elseif nen == 6  %P2
        N    = [xi.*(2*xi-1),et.*(2*et-1),(1-2*(xi+et)).*(1-(xi+et)),4*xi.*et,4*et.*(1-(xi+et)),4*xi.*(1-(xi+et))]; 
        Nxi  = [4*xi-1,zeros(size(xi)),-3+4*(xi+et),4*et,-4*et,4*(1-2*xi-et)]; 
        Neta = [zeros(size(xi)),4*et-1,-3+4*(xi+et),4*xi,4*(1-xi-2*et),-4*xi]; 
    elseif nen == 4  %P1+ 
        N    = [xi,et,1-(xi+et),27*xi.*et.*(1-xi-et)]; 
        Nxi  = [ones(size(xi)),zeros(size(xi)),-ones(size(xi)),27*et.*(1-2*xi-et)]; 
        Neta = [zeros(size(xi)),ones(size(xi)),-ones(size(xi)),27*xi.*(1-2*et-xi)]; 
    else
       error ('Error in ShapeFunc: unavailable triangle')
   end
else 
    error ('Error in ShapeFunc')
end
end
