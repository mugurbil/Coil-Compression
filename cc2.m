function [ new ] = cc2( data, n, p )
%CC2 coil compression 2D 
%   new_data = cc2( data ) 
%   reduces number of coils to virtual 
%   takes in raw k-space data
%   uses singular value decomposition if not p  

%   Mehmet Ugurbil, University of Minnesota, July 2012

tic

if nargin<3
    p=0;
end

if nargin<2 || isempty(n)
    n=6;
end

z=size(data);

%   create the data matrix
X=permute(data,[3 1 2]);
X=reshape(X, z(3), 1, []);
X=sq(X);

if p
    
    [c u]=princomp(X, 'econ');
    
else
    
    [u s ~]=svd(X, 'econ');
    
end


A=u(:,1:n)';
%   create the new data
new=A*X;
new=reshape(new, n, z(1), z(2));
new=permute(new, [2 3 1]);

toc

return

