function [X, T] = rectangleMesh(elem, nen, x1, x2, y1, y2, nx, ny)
% [X, T] = rectangleMesh(elem, nen, x1, x2, y1, y2, nx, ny)
% Input:
%    elem: element type (0: cuadrilateral,   1: triangle)
%    nen: number of element nodes
%    x1, x2, y1, y2: coordinates for a domain [x1, x2, y1, y2]
%    nx,ny: number of elements in each direction
% Output:
%    X: nodal coordinates
%    T: connectivities
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

npx = nx + 1;
npy = ny + 1;

X = zeros(npx*npy, 2);

xs = linspace(x1, x2, npx)';
unos = ones(npx, 1);
yys = linspace(y1, y2, npy);
for i = 1 : npy
    ys = yys(i)*unos;
    posi = (i - 1)*npx + 1 : i*npx;
    X(posi, :) = [xs, ys];
end

if elem == 0 % Cuadrilateral
    if nen == 4 % Q1
        T = zeros(nx*ny, 4);
        for i = 1 : ny
            for j = 1 : nx
                ielem = (i - 1)*nx + j;
                inode = (i - 1)*npx + j;
                T(ielem, :) = [inode, inode + 1, inode + (nx + 2), inode + npx];
            end
        end
    elseif nen == 9 %Q2
        if (nx - 2*floor(nx / 2) ~=0) && (ny - 2*floor(ny / 2) ~=0)
            error('Number of nodes in X or Y direction is not odd')
        end
        for i = 1 : ny / 2
            for j = 1 : nx / 2
                ielem = (i - 1)*nx / 2 + j;
                inode = (i - 1)*2*npx + 2*(j - 1) + 1;
                T(ielem, :) = [inode, inode + 2, inode + 2*npx + 2, inode + 2*npx, inode + 1, inode + 2 + npx, inode + 2*npx + 1, inode + npx, inode + npx + 1];
            end
        end
    end
elseif elem == 1   % Triangle
    if nen == 3 % P1
        T = zeros(nx*ny, 3);
        for i = 1 : ny
            for j = 1 : nx
                ielem = 2*((i - 1)*nx + j) - 1;
                inode = (i - 1)*npx + j;
                T(ielem, :) = [inode, inode + 1, inode + npx];
                T(ielem + 1, :) = [inode + 1, inode + 1 + npx, inode + npx];
            end
        end
        % Modificamos los elementos de las esquinas SO y NE para evitar que tengan todos los nodos con velocidades prescritas
        T(1, :) = [1, npx + 2, npx + 1];
        T(2, :) = [1, 2, npx + 2];
        aux = size(T, 1);
        T(aux, :) = [npx*ny - 1, npx*npy, npx*npy - 1];
        T(aux - 1, :) = [npx*ny - 1, npx*ny, npx*npy];
    elseif nen == 4 % P1+   Obs: no se genera el nodo interior (funcion burbuja)
        T = zeros(nx*ny, 4);
        for i = 1 : ny
            for j = 1 : nx
                ielem = 2*((i - 1)*nx + j) - 1;
                inode = (i - 1)*npx + j;
                n_ad = npx*npy + 2*((i - 1)*nx + j) - 1;
                T(ielem, :) = [inode, inode + 1, inode + npx, n_ad];
                T(ielem + 1, :) = [inode + 1, inode + 1 + npx, inode + npx, n_ad + 1];
            end
        end
        % Modificamos los elementos de las esquinas SO y NE para evitar que tengan todos los nodos con velocidades prescritas
        aux = size(T, 1);
        T(1, :) = [1, npx + 2, npx + 1, npx*npy + 1];
        T(2, :) = [1, 2, npx + 2, npx*npy + 2];
        T(aux - 1, :) = [npx*ny - 1, npx*npy, npx*npy - 1, npx*npy + 2*nx*ny - 1];
        T(aux, :)   = [npx*ny - 1, npx*ny, npx*npy, npx*npy + 2*nx*ny];
    elseif nen == 6  % P2
        for i = 1 : ny / 2
            for j = 1 : nx / 2
                ielem = 2*((i - 1)*nx / 2 + j) - 1;
                inode = (i - 1)*2*npx + 2*(j - 1) + 1;
                T(ielem, :) = [inode, inode + 2, inode + 2*npx, inode + 1, inode + 1 + npx, inode + npx];
                T(ielem + 1, :) = [inode + 2, inode + 2 + 2*npx, inode + 2*npx, inode + 2 + npx, inode + 1 + 2*npx, inode + 1 + npx];
            end
        end
        % Modificamos los elementos de las esquinas SO y NE para evitar que tengan todos los nodos con velocidades prescritas
        T(1, :) = [1, 2*npx + 3, 2*npx + 1, npx + 2, 2*npx + 2, npx + 1];
        T(2, :) = [1, 3, 2*npx + 3, 2, npx + 3, npx + 2];
        aux = size(T, 1);
        T(aux - 1, :) = [npx*(ny - 1) - 2, npx*(ny - 1), npx*npy, npx*(ny - 1) - 1, npx*ny, npx*ny - 1];
        T(aux, :) = [npx*(ny - 1) - 2, npx*npy, npx*npy - 2, npx*ny - 1, npx*npy - 1, npx*ny - 2 ];
    end
end