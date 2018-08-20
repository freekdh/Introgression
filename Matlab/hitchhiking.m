function hitchhiking()
    %(1)AB, (2)Ab, (3)aB, (4)ab
    
    WA = 1.1;
    Wa = 0.9;
    r = 0;
    size = 5;
    
    birthdeathmat = birthdeathv3(size,WA,Wa,r);
    %spy(birthdeathmat)
    %get probability of absorbption
    [V,D,W] = eig(full(birthdeathmat));
    absorbmatrix = sparse(size^4,size^4);
    for i = 1:size^4
        if(D(i,i)==1)
            elem = find(V(:,i) == 1);
            absorbmatrix(elem,:) = W(:,i)*(1/W(elem,i));
        end
    end
    
    absorbmatrix(:,mat2elem(1,3,3,1,5))
    
    save('absorbmat.mat', 'absorbmatrix');
end

function absorbtionmatrix()
    for n = 2:size
        %get rescue / hitchhiking but conditioned on fixation.
        vec = absorbmatrix(:,mat2elem(1,n,size,1,size));
        FA = 0;
        for i = 2:size^4 % skip the first state (0,0,0,0)
            if(vec(i) ~= 0)
                [AB,Ab,aB,ab] = elem2mat(i,size);
                [FAtemp] = projectF(AB-1,Ab-1,aB-1,ab-1,s,r,10);
                FA = FA+vec(i)*FAtemp/(1-vec(1)); %FAB/(FAB+FAb)
                %Fa = Fa+vec(i)*Fatemp/(1-vec(1));
            end
        end

        %FA
    end
end

function realizations(birthdeathmat)
        
    %realization
    mc = dtmc(transpose(birthdeathmat));
    t=100;
    out=simulate(mc,t,'X0',initialvec(1,3,size,1,size));
    matout = zeros(t,4);
    for i = 1:t
        [AB,Ab,aB,ab]=elem2mat(out(i),size);
        matout(i,:) = [AB,Ab,aB,ab];
    end
    
    plot(matout)
    dlmwrite("data/test.csv",matout)   
end

function FA = projectFAfromMatrix(size,mat)
    %get probability of absorbption
    [V,D,W] = eig(full(mat));
    absorbmatrix = sparse(size^4,size^4);
    for i = 1:size^4
        if(D(i,i)==1)
            elem = find(V(:,i) == 1);
            absorbmatrix(elem,:) = W(:,i)*(1/W(elem,i));
        end
    end
    
    %get rescue / hitchhiking but conditioned on fixation.
    vec = absorbmatrix(:,mat2elem(1,2,size,1,size));
    FA = 0;
    for i = 2:size^4 % skip the first state (0,0,0,0)
        if(vec(i) ~= 0)
            [AB,Ab,aB,ab] = elem2mat(i,size);
            [FAtemp] = projectF(AB-1,Ab-1,aB-1,ab-1,s,r,10);
            FA = FA+vec(i)*FAtemp/(1-vec(1)); %FAB/(FAB+FAb)
            %Fa = Fa+vec(i)*Fatemp/(1-vec(1));
        end
    end
    
    FA

end

function [E, V, MA, Ma] = data(t,initvec,powermat,size)
    resol = 100;
    
    temppower = eye(size^4);
    E = zeros(t,4);
    V = zeros(t,4);
    MA = zeros(t,resol);
    Ma = zeros(t,resol);
    
    for tt = 1:t
        vec = temppower*initvec;
        for i = 1: size^4
            [k,l,m,n] = elem2mat(i,size);
            E(tt,1) = E(tt,1) + (k-1)*vec(i);
            E(tt,2) = E(tt,2) + (l-1)*vec(i);
            E(tt,3) = E(tt,3) + (m-1)*vec(i);
            E(tt,4) = E(tt,4) + (n-1)*vec(i);
            V(tt,1) = V(tt,1) + (k-1)*(k-1)*vec(i);
            V(tt,2) = V(tt,2) + (l-1)*(l-1)*vec(i);
            V(tt,3) = V(tt,3) + (m-1)*(m-1)*vec(i);
            V(tt,4) = V(tt,4) + (n-1)*(n-1)*vec(i);
            
            %full prob dist for F
            if((k-1)+(l-1) ~= 0)
                %fixation of A occurred
                if(k-1 == 0)
                    MA(tt,1) = MA(tt,1) + vec(i);
                else
                    MA(tt,ceil(((k-1)/((k-1)+(l-1)))*resol)) =  MA(tt,ceil(((k-1)/((k-1)+(l-1)))*resol)) + vec(i);
                end
            else
                %no fixation of A occurred
                MA(tt,1) = MA(tt,1) + vec(i);
            end

            if((m-1)+(n-1) ~= 0)
                %fixation of A occurred
                if(m-1 == 0)
                    Ma(tt,1) = Ma(tt,1) + vec(i);
                else
                    Ma(tt,ceil(((m-1)/((m-1)+(n-1)))*resol)) =  Ma(tt,ceil(((m-1)/((m-1)+(n-1)))*resol)) + vec(i);
                end
            else
                %no fixation of A occurred
                Ma(tt,1) = Ma(tt,1) + vec(i);
            end
        end
        
        V(tt,1) = V(tt,1) - (E(tt,1)*E(tt,1));
        V(tt,2) = V(tt,2) - (E(tt,2)*E(tt,2));
        V(tt,3) = V(tt,3) - (E(tt,3)*E(tt,3));
        V(tt,4) = V(tt,4) - (E(tt,4)*E(tt,4));
                
        temppower = temppower*powermat;
    end
end

function vec = initialvec(AB,Ab,aB,ab,size)
vec = zeros(size^4,1);
vec(mat2elem(AB,Ab,aB,ab,size)) = 1;
end

function matrix = birthdeathv2(size, WA, Wa, r)
    matrix = sparse(size^4,size^4);
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    if(AB~=size && Ab~= size) %process ends if AB || Ab == size
                        sum=AB+Ab+aB+ab-4;
                        temp = 0;
                        if(sum~=0)
                            fAB=(AB-1)/sum;
                            fAb=(Ab-1)/sum;
                            faB=(aB-1)/sum;
                            fab=(ab-1)/sum;
                            D=(fAB*fab-fAb*faB);
                            dAB=fAB-r*D;
                            dAb=fAb+r*D;
                            daB=faB+r*D;
                            dab=fab-r*D;
                            %Calculate temp
                            if(AB ~= 1)
                            temp = temp + (AB-1)*WA*0.1;
                            end
                            if(Ab ~= 1)
                            temp = temp + (Ab-1)*WA*0.1;
                            end
                            if(aB ~= 1)
                            temp = temp + (aB-1)*Wa*0.1;
                            end
                            if(ab ~= 1)
                            temp = temp + (ab-1)*Wa*0.1;
                            end
                            if(AB ~= size)
                            temp = temp + (fAB-r*(fAB*fab-fAb*faB));
                            end
                            if(Ab ~= size)
                            temp = temp + (fAb+r*(fAB*fab-fAb*faB));
                            end
                            if(aB ~= size)
                            temp = temp + (faB+r*(fAB*fab-fAb*faB));
                            end
                            if(ab ~= size)
                            temp = temp + (fab-r*(fAB*fab-fAb*faB));
                            end

                            %actual probabilities conditioned on either birth
                            %or death happening.
                            if(temp ~= 0)                           
                                if(AB ~= 1)
                                matrix(mat2elem(AB-1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = fAB ; 
                                end
                                if(Ab ~= 1)
                                matrix(mat2elem(AB,Ab-1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (Ab-1)*WA*0.1/temp; 
                                end
                                if(aB ~= 1)
                                matrix(mat2elem(AB,Ab,aB-1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (aB-1)*Wa*0.1/temp; 
                                end
                                if(ab ~= 1)
                                matrix(mat2elem(AB,Ab,aB,ab-1,size),mat2elem(AB,Ab,aB,ab,size)) = (ab-1)*Wa*0.1/temp;
                                end

                                if(AB ~= size)
                                matrix(mat2elem(AB+1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((fAB-r*(fAB*fab-fAb*faB))/temp); 
                                end
                                if(Ab ~= size)
                                matrix(mat2elem(AB,Ab+1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((fAb+r*(fAB*fab-fAb*faB))/temp); 
                                end
                                if(aB ~= size)
                                matrix(mat2elem(AB,Ab,aB+1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((faB+r*(fAB*fab-fAb*faB))/temp); 
                                end
                                if(ab ~= size)
                                matrix(mat2elem(AB,Ab,aB,ab+1,size),mat2elem(AB,Ab,aB,ab,size)) = ((fab-r*(fAB*fab-fAb*faB))/temp); 
                                end
                            else
                                matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1;
                            end
                        end
                    else
                    %if AB == size || Ab == size
                    matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1;
                    end                    
                end
            end
        end
    end
    
    matrix(1,1) = 1;
    
end

function matrix = birthdeathv3(size, WA, Wa, r)
    matrix = sparse(size^4,size^4);
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    if(AB~=size && Ab~= size) %process ends if AB || Ab == size
                        sum=AB+Ab+aB+ab-4;
                        if(sum~=0)
                            fAB=(AB-1)/sum;
                            fAb=(Ab-1)/sum;
                            faB=(aB-1)/sum;
                            fab=(ab-1)/sum;
                            D=(fAB*fab-fAb*faB);
                            dAB=fAB-r*D;
                            dAb=fAb+r*D;
                            daB=faB+r*D;
                            dab=fab-r*D;
                            deathAB=2-WA;
                            deathAb=2-WA;
                            deathaB=2-Wa;
                            deathab=2-Wa;
                            weightedsumdeath = (fAB*deathAB+fAb*deathAb+faB*deathaB+fab*deathab);
                            Pbirth=1/(1+weightedsumdeath);
                            Pdeath=1-Pbirth;
                            if(AB ~= 1)
                            matrix(mat2elem(AB-1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pdeath*(fAB*deathAB)/(weightedsumdeath);
                            end
                            if(Ab ~= 1)
                            matrix(mat2elem(AB,Ab-1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pdeath*(fAb*deathAb)/(weightedsumdeath); 
                            end
                            if(aB ~= 1)
                            matrix(mat2elem(AB,Ab,aB-1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pdeath*(faB*deathaB)/(weightedsumdeath);
                            end
                            if(ab ~= 1)
                            matrix(mat2elem(AB,Ab,aB,ab-1,size),mat2elem(AB,Ab,aB,ab,size)) = Pdeath*(fab*deathab)/(weightedsumdeath);
                            end

                            matrix(mat2elem(AB+1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pbirth*dAB; 
                            matrix(mat2elem(AB,Ab+1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pbirth*dAb;                                 
                            if(aB ~= size)
                                matrix(mat2elem(AB,Ab,aB+1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pbirth*daB; 
                            else
                                matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = Pbirth*daB;                                 
                            end
                            if(ab ~= size)
                                matrix(mat2elem(AB,Ab,aB,ab+1,size),mat2elem(AB,Ab,aB,ab,size)) = Pbirth*dab; 
                            else
                                matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) + Pbirth*dab;                                 
                            end
                        else
                            matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1;
                        end
                    else
                    %if AB == size || Ab == size
                    matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1;
                    end                    
                end
            end
        end
    end
    matrix(1,1) = 1;
    
end

function index = mat2elem(i,j,k,l,N)
index = (l-1)*N*N*N+(k-1)*N*N+(j-1)*N+i;
end

function [i,j,k,l] = elem2mat(index, N)
    i = mod(index-1,N)+1;
    j = mod(floor((index-1)/(N)),N)+1;
    k = mod(floor((index-1)/(N*N)),N)+1;
    l = mod(floor((index-1)/(N*N*N)),N)+1;
end

function [FA] = projectF(AB,Ab,aB,ab,s,r,t)
    sum = AB+Ab+aB+ab;
    fAB = AB/sum;
    fAb = Ab/sum;
    faB = aB/sum;
    fab = ab/sum;
    for i = 1:t
        Wbar = fAB*(1+s)+fAb*(1+s)+faB+fab;
        fABn=(fAB*fAB*(1+s)^2+fAB*faB*(1+s)+fAB*fAb*(1+s)^2+fAB*fab*(1+s)*(1-r)+faB*fAb*(1+s)*r)/(Wbar^2);
        fAbn=(fAb*fAb*(1+s)^2+fAb*fab*(1+s)+fAb*fAB*(1+s)^2+fAB*fab*(1+s)*r+faB*fAb*(1+s)*(1-r))/(Wbar^2);
        faBn=(faB*faB+faB*fab+faB*fAB*(1+s)+fAB*fab*(1+s)*r+faB*fAb*(1+s)*(1-r))/(Wbar^2);
        fabn=(fab*fab+fab*faB+fab*fAb*(1+s)+fAB*fab*(1+s)*(1-r)+faB*fAb*(1+s)*r)/(Wbar^2);
        fAB = fABn;
        fAb = fAbn;
        faB = faBn;
        fab = fabn;
    end
    FA = fAB/(fAB+fAb);
    %Fa = faB/(faB+fab);
end

% CODE GRAVEYARD

function [i,j,k,l] = multinomial(AB,Ab,aB,ab,ABn,Abn,aBn,abn,wA,wa,r)
    sum = AB+Ab+aB+ab;
    if(sum ~= 0)
        %initialize
        fAB=AB/sum;
        fAb=Ab/sum;
        faB=aB/sum;
        fab=ab/sum;
        Wbar=(fAB+fAb)*wA+(faB+fab)*wa;

        %selection
        fABs = fAB*wA/Wbar;
        fAbs = fAb*wA/Wbar;
        faBs = faB*wa/Wbar;
        fabs = fab*wa/Wbar;

        %recombination
        D = fABs*fabs-fAbs*faBs;
        fABr = fABs - r*D;
        fAbr = fAbs + r*D;
        faBr = faBs + r*D;
        fabr = fabs - r*D;

        %multinomial
        fABr^ABn*fAbr^Abn*faB^aBn*faBr^aBn*fabr^abn*nchoosek(n,k);

        
    end
    
end

function matrix = recombinationMatrix(size, rec,dt)
    matrix = sparse(size^4,size^4);
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    sum=AB+Ab+aB+ab-4;
                    temp =0;
                    if(AB ~= size && Ab ~= 1 && aB ~= 1 && ab ~= size)
                    matrix(mat2elem(AB+1,Ab-1,aB-1,ab+1,size), mat2elem(AB,Ab,aB,ab,size)) = (Ab-1)*(aB-1)*dt*rec; temp = temp + (Ab-1)*(aB-1)*dt*rec;
                    end
                    if(AB ~= 1 && Ab ~= size && aB ~= size && ab ~= 1)
                    matrix(mat2elem(AB-1,Ab+1,aB+1,ab-1,size), mat2elem(AB,Ab,aB,ab,size)) = (AB-1)*(ab-1)*dt*rec; temp = temp + (AB-1)*(ab-1)*dt*rec;                 
                    end
                    
                    if(sum ~= 0)
                    matrix(mat2elem(AB,Ab,aB,ab,size), mat2elem(AB,Ab,aB,ab,size)) = 1-temp;
                    end
                end
            end
        end
    end
    
    matrix(1,1) = 1;
    
end

function matrix = selectionMatrix(size, bA, dA, ba, da, dt)
    matrix = sparse(size^4,size^4);
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    sum=AB+Ab+aB+ab-4;
                    temp = 0;
                    if(AB ~= 1)
                    matrix(mat2elem(AB-1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (AB-1)*dt*dA; temp=temp+(AB-1)*dt*dA;
                    end
                    if(Ab ~= 1)
                    matrix(mat2elem(AB,Ab-1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (Ab-1)*dt*dA; temp=temp+(Ab-1)*dt*dA;
                    end
                    if(aB ~= 1)
                    matrix(mat2elem(AB,Ab,aB-1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (aB-1)*dt*da; temp=temp+(aB-1)*dt*da;
                    end
                    if(ab ~= 1)
                    matrix(mat2elem(AB,Ab,aB,ab-1,size),mat2elem(AB,Ab,aB,ab,size)) = (ab-1)*dt*da; temp=temp+(ab-1)*dt*da;
                    end
                    
                    if(AB ~= 1 && AB ~= size)
                    matrix(mat2elem(AB+1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (AB-1)*dt*bA; temp=temp+(AB-1)*dt*bA;
                    end
                    if(Ab ~= 1 && Ab ~= size)
                    matrix(mat2elem(AB,Ab+1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (Ab-1)*dt*bA; temp=temp+(Ab-1)*dt*bA;
                    end
                    if(aB ~= 1 && aB ~= size)
                    matrix(mat2elem(AB,Ab,aB+1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (aB-1)*dt*ba; temp=temp+(aB-1)*dt*ba;
                    end
                    if(ab ~= 1 && ab ~= size)
                    matrix(mat2elem(AB,Ab,aB,ab+1,size),mat2elem(AB,Ab,aB,ab,size)) = (ab-1)*dt*ba; temp=temp+(ab-1)*dt*ba;
                    end
                    
                    
                    if(sum ~= 0)
                    matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1-temp;
                    end
                end
            end
        end
    end
    
    matrix(1,1) = 1;
    
end

function matrix = birthdeathv1(size, sA, r)
    matrix = sparse(size^4,size^4);
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    sum=AB+Ab+aB+ab-4;
                    temp = 0;
                    if(sum~=0)
                        fAB=(AB-1)/sum;
                        fAb=(Ab-1)/sum;
                        faB=(aB-1)/sum;
                        fab=(ab-1)/sum;
                        %Calculate temp
                        if(AB ~= 1 && AB ~= size)
                        temp = temp + (AB-1)*(1-sA)*0.1;
                        end
                        if(Ab ~= 1 && Ab ~= size)
                        temp = temp + (Ab-1)*(1-sA)*0.1;
                        end
                        if(aB ~= 1)
                        temp = temp + (aB-1)*0.1;
                        end
                        if(ab ~= 1)
                        temp = temp + (ab-1)*0.1;
                        end
                        if(AB ~= size)
                        temp = temp + (fAB-r*(fAB*fab-fAb*faB));
                        end
                        if(Ab ~= size)
                        temp = temp + (fAb+r*(fAB*fab-fAb*faB));
                        end
                        if(aB ~= size)
                        temp = temp + (faB+r*(fAB*fab-fAb*faB));
                        end
                        if(ab ~= size)
                        temp = temp + (fab-r*(fAB*fab-fAb*faB));
                        end
                        
                        %actual probabilities conditioned on either birth
                        %or death happening.
                        if(temp ~= 0)
                            if(AB ~= 1 && AB ~= size)
                            matrix(mat2elem(AB-1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (AB-1)*(1-sA)*0.1/temp; 
                            end
                            if(Ab ~= 1 && Ab ~= size)
                            matrix(mat2elem(AB,Ab-1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (Ab-1)*(1-sA)*0.1/temp; 
                            end
                            if(aB ~= 1)
                            matrix(mat2elem(AB,Ab,aB-1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = (aB-1)*0.1/temp; 
                            end
                            if(ab ~= 1)
                            matrix(mat2elem(AB,Ab,aB,ab-1,size),mat2elem(AB,Ab,aB,ab,size)) = (ab-1)*0.1/temp;
                            end

                            if(AB ~= size)
                            matrix(mat2elem(AB+1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((fAB-r*(fAB*fab-fAb*faB))/temp); 
                            end
                            if(Ab ~= size)
                            matrix(mat2elem(AB,Ab+1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((fAb+r*(fAB*fab-fAb*faB))/temp); 
                            end
                            if(aB ~= size)
                            matrix(mat2elem(AB,Ab,aB+1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((faB+r*(fAB*fab-fAb*faB))/temp); 
                            end
                            if(ab ~= size)
                            matrix(mat2elem(AB,Ab,aB,ab+1,size),mat2elem(AB,Ab,aB,ab,size)) = ((fab-r*(fAB*fab-fAb*faB))/temp); 
                            end
                        else
                            matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1;
                        end
                    end
                end
            end
        end
    end
    
    matrix(1,1) = 1;
    
end


