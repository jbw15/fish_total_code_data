function [] = func_inputfiles_creation(points)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
 [NNum nnod a]= size(points);  %[718 1782];  
 
 %NNum = 6000; 
 start = 1;

for n = start:NNum
     n1 = n-start+1;
    if (n1<10)
         fn = ['0/s.000000' num2str(n-start+1)];
    elseif(n1<100)
         fn = ['0/s.00000' num2str(n-start+1)];
    elseif(n1<1000)
         fn = ['0/s.0000' num2str(n-start+1)];
    elseif(n1<10000)
         fn= ['0/s.000' num2str(n-start+1)];
    else
        fn= ['0/s.00' num2str(n-start+1)];
    end
    
    ss=zeros(nnod,3);
         for i=1:nnod
             for j=1:3
             ss(i,j)=points(n,i,j);
             %ss(i,j+3)=vels(n,i,j);
             end
         end
         
         dlmwrite(fn,ss,'delimiter','\t');
end
 
end





