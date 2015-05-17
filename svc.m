function [ PA ] = svc( A )
%SVC smooth virtual coils
%   add unitary P matrix
%   get smoother virtual coil sensitivities 
%   newA=svc(A);

%   Mehmet Ugurbil, University of Minnesota, July 2012

disp('Virtual Coil Alignment In Progress');
t=tic;

%   intialize the matrix to be returned
PA=zeros(size(A));

%   set P1=I or A1_new=A1_old
PA(1,:,:)=sq(A(1,:,:));

%   calculate and apply P matrices for remaining indecies
for x=2:size(A,1)
    
    %   define nvc*nvc matrix C
    C=sq(A(x,:,:))*sq(PA(x-1,:,:))';
   
     %   single value decomposition
    [u s v]=svd(C,'econ');

    P=v*u';
    
    PA(x,:,:)=P*sq(A(x,:,:));
    
    
end

toc(t)
return