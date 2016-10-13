function Be = Bmat2D(nnel, dNdx, dNdy)
% Be = Bmat2D(nnel, dNdx, dNdy)
%------------------------------------------------------------------------
%  Purpose:
%     determine the strain-displacement matrix for 2D solids
%
%  Synopsis:
%     Be = get_Bmat_2D(nnel,dNdx,dNdy)
%
%  Variable Description:
%       Be  : strain-displacement matrix of element
%     nnel  : number of nodes per element
%     dNdx  : derivatives of shape functions with respect to x
%     dNdy  : derivatives of shape functions with respect to y
%------------------------------------------------------------------------
%      [dNdx(1) 0       dNdx(2) 0       ... dNdx(nnel) 0         ]
% Be = [0       dNdy(1) 0       dNdy(2) ... 0          dNdy(nnel)]
%      [dNdy(1) dNdx(1) dNdy(2) dNdx(2) ... dNdy(nnel) dNdx(nnel)]
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

Be = zeros(3, 2*nnel);
for i = 1 : nnel
    i1 = (i - 1)*2 + 1;
    i2 = i1 + 1;
    Be(1, i1) = dNdx(i);
    Be(2, i2) = dNdy(i);
    Be(3, i1) = dNdy(i);
    Be(3, i2) = dNdx(i);
end
end
