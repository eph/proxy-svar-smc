function [nbplt,nr,nc,lr,lc,nstar] = pltorg(number)
% stephane.adjemian@cepremap.cnrs.fr [06-07-2004]
nrstar = 3;
ncstar = 2;
nstar  = nrstar*ncstar;
nbplt  = 0;
nr     = 0;
nc     = 0;
lr     = 0;
lc     = 0;
if number == 1
    nbplt = 1;
    nr    = 1;
    nc    = 1;
elseif number == 2
    nbplt = 1;
    nr    = 2;
    nc    = 1;
elseif number == 3
    nbplt = 1;
    nr    = 2;
    nc    = 2;
elseif number == 4
    nbplt = 1;
    nr    = 2;
    nc    = 2;
elseif number == 5
    nbplt = 1;
    nr    = 2;
    nc    = 3;
elseif number == 6
    nbplt = 1;
    nr    = 2;
    nc    = 3;
elseif number == 8
    nbplt = 1;
    nr    = 2;
    nc    = 4;
elseif number == 7
    nbplt = 1;
    nr    = 4;
    nc    = 2;
elseif number == 9
    nbplt = 1;
    nr    = 3;
    nc    = 3;
elseif number == 10
    nbplt = 1;
    nr    = 5;
    nc    = 2;
elseif number == 12
    nbplt = 1;
    nr    = 4;
    nc    = 3;
elseif number == 13
    nbplt = 1;
    nr    = 4;
    nc    = 4;
elseif number == 14
    nbplt = 1;
    nr    = 4;
    nc    = 4;
else
    if number/nstar == round(number/nstar)
        nbplt = number/nstar;
        nr    = nrstar;
        nc    = ncstar;
        lr    = nr;
        lc    = nc; 
    else
        nbplt = ceil(number/nstar);
        nr    = nrstar;
        nc    = ncstar;
        reste = number-(nbplt-1)*nstar;
        if reste == 1
            lr    = 1;
            lc    = 1;
        elseif reste == 2
            lr    = 2;
            lc    = 1;
        elseif reste == 3
            lr    = 3;
            lc    = 1;
        elseif reste == 4
            lr    = 2;
            lc    = 2;
        elseif reste == 5
            lr    = 3;
            lc    = 2;
        elseif reste == 6
            lr    = 3;
            lc    = 2;    
        elseif reste == 7
            lr    = 3;
            lc    = 3;    
        elseif reste == 8
            lr    = 3;
            lc    = 3;
        end
    end
end
