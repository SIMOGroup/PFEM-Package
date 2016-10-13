function [pospg, wpg] = Quadrature(elem, degreeGauss)
% [pospg, wpg] = Quadrature(elem, degreeGauss)
%
% Input:
%   elem:   type of element (0 for quadrilaterals and 1 for triangles)
%   ngaus:  number of Gauss points in each element
% Output:
%   pospg, wpg: Gauss points and weights on the reference element
%
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

if elem == 0
    pos1 = 1/sqrt(3);
    pospg=[-pos1   -pos1
        pos1   -pos1
        pos1    pos1
        -pos1    pos1];
    wpg=[ 1 1 1 1];
elseif elem == 1
    if degreeGauss == 1  % degree 1
        pos1 = 1/3;
        pospg=[pos1   pos1];
        wpg = 1/2;
    elseif degreeGauss == 2  % degree 2
        pospg = [1/6   1/6
            4/6   1/6
            1/6   4/6];
        pes1 = 1/6;
        wpg = [pes1   pes1   pes1];
    elseif degreeGauss == 3  % degree 3
        pospg=[1/3   1/3
            0.6   0.2
            0.2   0.6
            0.2   0.2];
        wpg = [-27/96   25/96   25/96   25/96];
    elseif degreeGauss == 5   % degree 5
        a = 0.101286507323456338800987361915123;
        b = 0.470142064105115089770441209513447;
        P1 = 0.0629695902724135762978419727500906;
        P2 = 0.0661970763942530903688246939165759;
        pospg=[a       a
            a       1-2*a
            1-2*a   a
            b       b
            b       1-2*b
            1-2*b   b
            1/3     1/3];
        wpg = [P1, P1, P1, P2, P2, P2, 0.1125];
    elseif degreeGauss == 6   % degree 6
        a = 0.0630890144915022283403316028708191;
        b = 0.249286745170910421291638553107019;
        c = 0.0531450498448169473532496716313981;
        d = 0.310352451033784405416607733956552;
        P1 = 0.0254224531851034084604684045534344;
        P2 = 0.0583931378631896830126448056927897;
        P3 = 0.0414255378091867875967767282102212;
        pospg = [a        a
            a        1-2*a
            1-2*a    a
            b        b
            b        1-2*b
            1-2*b    b
            c        d
            c        1-c-d
            1-c-d    c
            1-c-d    d
            d        1-c-d
            d        c  ];
        wpg = [P1, P1, P1, P2, P2, P2, P3, P3, P3, P3, P3, P3];
    elseif degreeGauss == 8  % degree 8
        a = 0.170569307751760206622293501491464;
        b = 0.0505472283170309754584235505965989;
        c = 0.459292588292723156028815514494169;
        d = 0.728492392955404281241000379176061;
        e = 0.263112829634638113421785786284643;
        f = 1/3;
        P1 = 0.0516086852673591251408957751460645;
        P2 = 0.0162292488115990401554629641708902;
        P3 = 0.0475458171336423123969480521942921;
        P4 = 0.0136151570872174971324223450369544;
        P5 = 0.0721578038388935841255455552445323;
        pospg = [a        a
            a        1-2*a
            1-2*a    a
            b        b
            b        1-2*b
            1-2*b    b
            c        c
            c        1-2*c
            1-2*c    c
            d        e
            d        1-d-e
            1-d-e    d
            1-d-e    e
            e        1-d-e
            e        d
            f        f ];
        wpg = [P1*ones(1,3), P2*ones(1,3), P3*ones(1,3), P4*ones(1,6), P5];
    else
        error('Quadrature no contemplada')
    end
else
    error(' Quadrature no contemplada')
end
end
