function hitchhiking1D()
    
    WA = 1.01;
    size = 1000;
    
    mat = birthdeath(size,WA);
    
    initvec=initialvec(2,size);
    
    E = data(200,initvec,mat,size);
    save('growth1D.mat', 'E');
end

function vec = initialvec(AB,size)
vec = zeros(size,1);
vec(AB) = 1;
end

function matrix = birthdeath(size, W)
    matrix = sparse(size,size);
    for AB = 1:size
        if(AB~=size) %process ends if AB || Ab == size
            pcd = 0.3;
            birthrate = (AB-1)*((W-1)+pcd);
            deathrate = (AB-1)*pcd;
            if(AB ~= 1)
            matrix(AB-1,AB) = deathrate/(deathrate+birthrate);
            matrix(AB+1,AB) = birthrate/(deathrate+birthrate); 
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
        for AB = 2: size
            E(tt) = E(tt) + (AB-1)*(vec(AB)/(1-vec(1)));
        end        
        temppower = temppower*powermat;
    end
    
end
