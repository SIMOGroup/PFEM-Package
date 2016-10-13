function index = assembly(nd, nnel, ndof)
%----------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element
%
%  Synopsis:
%     [index]=feeldof(nd, nnel, ndof)
%
%  Variable Description:
%     index - system dof vector associated with element "iel"
%     iel - element number whose system dofs are to be determined
%     nnel - number of nodes per element
%     ndof - number of dofs per node
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

index = zeros(1, nnel*ndof);
k = 0;
for i = 1 : nnel
    start = (nd(i) - 1)*ndof;
    for j = 1 : ndof
        k = k + 1;
        index(k) = start + j;
    end
end

end

