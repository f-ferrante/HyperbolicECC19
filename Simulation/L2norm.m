function [ output_args ]=SpatialNorm(input_args)
Y=zeros(1,1);
for i=1:max(size(V1to))
    
   for j=1:min(size(V1to))
        
   Y(j)=((V1to(i,j).^2+V2to(i,j).^2)); 
    
    end
    Int(i)=sqrt(trapz(x,Y));
end
end