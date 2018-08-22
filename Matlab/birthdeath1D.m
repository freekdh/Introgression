function hitchhiking1D()
    
    WA = 1.01;
    size = 100;
    
    mat = birthdeath(size,WA,0.01);
        
    initvec=initialvec(2,size);
    E = data(1000,initvec,mat,size);
    E;
    semilogy(E)
end

function vec = initialvec(AB,size)
vec = zeros(size,1);
vec(AB) = 1;
end

function matrix = birthdeath(size, W,dt)
    matrix = sparse(size,size);
    for AB = 1:size
        if(AB~=size) %process ends if AB || Ab == size
            pcd = 0.3;
            birthrate = (AB-1)*((W-1)+pcd)*dt;
            deathrate = (AB-1)*pcd*dt;
            if(AB ~= 1)
            matrix(AB-1,AB) = deathrate;
            matrix(AB,AB)   = 1-deathrate-birthrate; 
            matrix(AB+1,AB) = birthrate; 
            else
            matrix(AB,AB) = 1;
            end
        else
        %if AB == size || Ab == size
        matrix(AB,AB) = 1;
        end                    
    end
end

function E = data(t,initvec,powermat,size)
    temppower = eye(size);
    E = zeros(t,1);
    
    for tt = 1:t
        vec = temppower*initvec;
        temp = 0;
        for AB = 2: size
            temp = temp + (AB-1)*(vec(AB)/(1-vec(1)));
        end    
        E(tt,1)=temp;
        temppower = temppower*powermat;
    end
    
end
