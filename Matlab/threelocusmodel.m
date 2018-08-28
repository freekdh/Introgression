function threelocusmodel ()
    
    p0 = 0.1;
    initvec = [p0,p0,p0,p0*(1-p0),p0*(1-p0),p0*(1-p0),p0-3*p0^2+2*p0^3];
    
    % All dynamics:
    %out = projectall(50, initvec,0.2,-0.1,0.5,0.5,0.5);
    
    nsk=-(1:10)/100;
    nsi=(1:10)/10;
    outmat = zeros(10,2);
    for sk=1:100
        out = projectend(50, initvec,0.3,-sk/1000,0.5,0.5,0.5);
        outmat(sk,1) = -sk/1000;
        outmat(sk,2) = projectFAend(out);
    end
    
    outmat
end

function [outlist] = itterate(inlist,si,sk,rij,rjk,rik)
    pi = inlist(1);
    pj = inlist(2);
    pk = inlist(3);
    dij = inlist(4);
    djk = inlist(5);
    dik = inlist(6);
    dijk = inlist(7);
    
    wbar = 1+si*pi+sk*pk;
    
    npi = (si*pi*(1-pi)+sk*dik)/(wbar)+pi;
    npj = (si*dij+sk*djk)/(wbar)+pj;
    npk = (sk*pk*(1-pk)+si*dik)/(wbar)+pk;
    ndij = (si*(1-2*pi)*dij+sk*dijk)/wbar + dij;
    ndjk = (sk*(1-2*pk)*djk+si*dijk)/wbar + djk;
    ndik = (si*((1-2*pi)*dik)+sk*((1-2*pk)*dik))/wbar + dik;
    ndijk = (si*(pi*(1-pi)*djk+(1-2*pi)*dijk)+sk*(pk*(1-pk)*dij+(1-2*pk)*dijk))/wbar + dijk;
    
    nnpi = npi;
    nnpj = npj;
    nnpk = npk;
    nndij = ndij-(pi-npi)*(pj-npj);
    nndjk = ndjk-(pj-npj)*(pk-npk);
    nndik = ndik-(pi-npi)*(pk-npk);
    nndijk = ndijk + ndij*(pk-npk)+ndik*(pj-npj)+ndjk*(pi-npi)-(pi-npi)*(pj-npj)*(pk-npk);
    
    outlist = zeros(1,7);

    outlist(1,1) = nnpi;
    outlist(1,2) = nnpj;
    outlist(1,3) = nnpk;
    outlist(1,4) = nndij*(1-rij);
    outlist(1,5) = nndjk*(1-rjk);
    outlist(1,6) = nndik*(1-rik);
    outlist(1,7) = nndijk*(1-rij)*(1-rjk);
    
end

function [outlist] = projectall(t, initlist, si, sk, rij, rjk, rik)
    
    outlist = zeros(t,8);
    templist = initlist;
    outlist(1,:) = genotypes(templist);
    for i = 2:t
        templist = itterate(templist,si,sk,rij,rjk,rik);
        outlist(i,:) = genotypes(templist);
    end
end

function [outlist] = projectend(t, initlist, si, sk ,rij, rjk, rik)
 
    outlist = zeros(1,8);
    templist = initlist;
    for i = 2:t
        templist = itterate(templist,si,sk,rij,rjk,rik);
    end
    outlist(1,:) = genotypes(templist);
end

function [outlist] = genotypes(inlist)
outlist = zeros(1,8);
outlist(1) = genotype(inlist,[0,0,0]);
outlist(2) = genotype(inlist,[0,0,1]);
outlist(3) = genotype(inlist,[0,1,0]);
outlist(4) = genotype(inlist,[0,1,1]);
outlist(5) = genotype(inlist,[1,0,0]);
outlist(6) = genotype(inlist,[1,0,1]);
outlist(7) = genotype(inlist,[1,1,0]);
outlist(8) = genotype(inlist,[1,1,1]);

end

function f = genotype(metric, X)
    f=(X(1)*metric(1)+(1-X(1))*(1-metric(1)))* ...
    (X(2)*metric(2)+(1-X(2))*(1-metric(2)))* ...
    (X(3)*metric(3)+(1-X(3))*(1-metric(3)))+ ...
    metric(4)*(-1)^(2-X(1)-X(2))*(X(3)*metric(3)+(1-X(3))*(1-metric(3))) + ...
    metric(5)*(-1)^(2-X(2)-X(3))*(X(1)*metric(1)+(1-X(1))*(1-metric(1))) + ...
    metric(6)*(-1)^(2-X(1)-X(3))*(X(2)*metric(2)+(1-X(2))*(1-metric(2))) + ...
    metric(7)*(-1)^(3-X(1)-X(2)-X(3));
end

function [FA,Faa] = calcFA(inlist)
    FA = (inlist(5)+inlist(6)) / (inlist(5)+inlist(6)+inlist(7)+inlist(8));
    Faa = (inlist(1)+inlist(2)) / (inlist(1)+inlist(2)+inlist(3)+inlist(4));    
end

function [outlist] = projectFAall(genotypesoutput)
    m = size(genotypesoutput);
    outlist = zeros(m(1),2);
    for t = 1:m(1)
        [a,b]=calcFA(genotypesoutput(t,:));
        outlist(t,1) = a;
        outlist(t,2) = b;
    end
end

function [outlist] = projectFAend(genotypesoutput)
    [a,b]=calcFA(genotypesoutput);
    outlist = a;
end

