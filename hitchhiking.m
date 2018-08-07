function hitchhiking()
    
    size = 6;
    mat = recombinationMatrix(0.5,size); % don't go over 12
    mat2 = selectionMatrix(size); % don't go over 12
    
    initvec = initialvec(1,2,6,1,size);
    powermat = mat*mat2;
        
    [AB,Ab,aB,ab] = expected(2,initvec,powermat,size);
    [AB,Ab,aB,ab]
end

function [AB,Ab,aB,ab] = expected(t, initvec, powermat, size)

    vec = powermat^t*initvec;

    AB = 0;
    Ab = 0;
    aB = 0;
    ab = 0;
    
    for i = 1: size^4
        [k,l,m,n] = elem2mat(i,size);
        AB = AB + k*vec(i);
        Ab = Ab + l*vec(i);
        aB = aB + m*vec(i);
        ab = ab + n*vec(i);
    end
    
end

function vec = initialvec(AB,Ab,aB,ab,size)
vec = zeros(size^4,1);
vec(mat2elem(AB,Ab,aB,ab,size)) = 1;
end

function matrix = recombinationMatrix(rec,size)
    matrix = sparse(size^4,size^4);
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    sum=AB+Ab+aB+ab-4;
                    temp =0;
                    if(AB ~= size && Ab ~= 1 && aB ~= 1 && ab ~= size)
                    matrix(mat2elem(AB+1,Ab-1,aB-1,ab+1,size), mat2elem(AB,Ab,aB,ab,size)) = ((Ab-1)/sum)*((aB-1)/sum)*rec; temp = temp + ((Ab-1)/sum)*((aB-1)/sum)*rec;
                    end
                    if(AB ~= 1 && Ab ~= size && aB ~= size && ab ~= 1)
                    matrix(mat2elem(AB-1,Ab+1,aB+1,ab-1,size), mat2elem(AB,Ab,aB,ab,size)) = ((AB-1)/sum)*((ab-1)/sum)*rec; temp = temp + ((AB-1)/sum)*((ab-1)/sum)*rec;                 
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

function matrix = selectionMatrix(size)
    matrix = sparse(size^4,size^4);
    bA = 0.2;
    dA = 0.02;
    ba = 0.02;
    da = 0.2;
    for AB = 1:size
        for Ab = 1:size
            for aB = 1:size
                for ab = 1:size
                    sum=AB+Ab+aB+ab-4;
                    temp = 0;
                    if(AB ~= 1)
                    matrix(mat2elem(AB-1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((AB-1)/sum)*dA; temp=temp+((AB-1)/sum)*dA;
                    end
                    if(Ab ~= 1)
                    matrix(mat2elem(AB,Ab-1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((Ab-1)/sum)*dA; temp=temp+((Ab-1)/sum)*dA;
                    end
                    if(aB ~= 1)
                    matrix(mat2elem(AB,Ab,aB-1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((aB-1)/sum)*da; temp=temp+((aB-1)/sum)*da;
                    end
                    if(ab ~= 1)
                    matrix(mat2elem(AB,Ab,aB,ab-1,size),mat2elem(AB,Ab,aB,ab,size)) = ((ab-1)/sum)*da; temp=temp+((ab-1)/sum)*da;
                    end
                    
                    if(AB ~= 1 && AB ~= size)
                    matrix(mat2elem(AB+1,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((AB-1)/sum)*bA; temp=temp+((AB-1)/sum)*bA;
                    end
                    if(Ab ~= 1 && Ab ~= size)
                    matrix(mat2elem(AB,Ab+1,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((Ab-1)/sum)*bA; temp=temp+((Ab-1)/sum)*bA;
                    end
                    if(aB ~= 1 && aB ~= size)
                    matrix(mat2elem(AB,Ab,aB+1,ab,size),mat2elem(AB,Ab,aB,ab,size)) = ((aB-1)/sum)*ba; temp=temp+((aB-1)/sum)*ba;
                    end
                    if(ab ~= 1 && ab ~= size)
                    matrix(mat2elem(AB,Ab,aB,ab+1,size),mat2elem(AB,Ab,aB,ab,size)) = ((ab-1)/sum)*ba; temp=temp+((ab-1)/sum)*ba;
                    end
                    
                    
                    if(sum ~= 0)
                    matrix(mat2elem(AB,Ab,aB,ab,size),mat2elem(AB,Ab,aB,ab,size)) = 1-temp;
                    end
                end
            end
        end
    end
    
    matrix(1,1) = 1;
    
    % now the upper limit 
    
    
end

function index = mat2elem(i,j,k,l,N)
index = (l-1)*N*N*N+(k-1)*N*N+(j-1)*N+i;
end

function [i,j,k,l] = elem2mat (index, N)
    i = mod(index-1,N)+1;
    j = mod(floor((index-1)/(N)),N)+1;
    k = mod(floor((index-1)/(N*N)),N)+1;
    l = mod(floor((index-1)/(N*N*N)),N)+1;
end

function testfunc(size)
    matrix = sparse(size^4,size^4);
    matrix(mat2elem(2,2,2,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(3,2,2,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(2,3,2,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(2,2,3,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(2,2,2,3,size),mat2elem(2,2,2,2,size)) = 1;

    matrix(mat2elem(1,2,2,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(2,1,2,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(2,2,1,2,size),mat2elem(2,2,2,2,size)) = 1;
    matrix(mat2elem(2,2,2,1,size),mat2elem(2,2,2,2,size)) = 1; 
    
    matrix(mat2elem(3,3,3,3,size),mat2elem(2,2,2,2,size)) = 1;    
    matrix(mat2elem(1,1,1,1,size),mat2elem(2,2,2,2,size)) = 1;    
    
    spy(matrix) 
end

