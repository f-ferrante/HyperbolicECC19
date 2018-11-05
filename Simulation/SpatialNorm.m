%This function generates a vector whose entries correspond to the (spatial) L2 norm of the solution (V1to(t,.),
%V2to(t,.)) 

function Int=SpatialNorm(V1to, V2to, x)

Y=zeros(1,1);
for i=1:max(size(V1to))
    
   for j=1:min(size(V1to))
        
   Y(j)=((V1to(i,j).^2+V2to(i,j).^2)); 
    
    end
    Int(i)=sqrt(trapz(x,Y));
end
end