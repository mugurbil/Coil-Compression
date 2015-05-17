function [ new ] = cc( data,n,smooth,p )
%CC coil compression 
%   reduces number of coils by creating virtual coils
%   takes fft along readout direction (x)
%   data: raw k-space data
%   n: optional, number of virtual coils, default is 6
%   smooth: optional, 0 if you want to avoid smoothing
%   p: optional, 1 to use pca instead of svd

%   Mehmet Ugurbil, University of Minnesota, July 2012

tic;

%   initialize parameters
if nargin<4
    p=0;
end

if nargin<3
    smooth=1;
end

if nargin<2
    n=6;
end

z=size(data);

%   take fft along readout direction(x)
temp=fftshift(fft(fftshift(data)));

%   A holds all the Ax matrices
A=zeros(z(1),n,z(3));

disp('Geometric Decomposition In Progress');

%   X is the data matrix
X=permute(temp,[1 3 2]);

t=tic;
for x=1:z(1)
    
    % %  geometric decomposition for coil compression algorithm % %
    
    %   Xx is the data matrix for single x value
    %   its rows consist of all v(k_y) values from a single coil
    Xx=sq(X(x,:,:));
    if p
         %   principle componant analysis
        [c u]=princomp(Xx,'econ');
    else
         %   single value decomposition
        [u s ~]=svd(Xx,'econ');
    end
    
    %   A is the first nvc rows from U'
    %   size(A) = nvc*nc
    Atemp=u(:,1:n)';

    A(x,:,:)=Atemp;
    
    
end
toc(t)

if smooth
    
    % %  smoothing by virtual coil alignment algorithm % %
    A=svc(A);
    
end

%   this will be the returned data set
new=zeros(z(1),n,z(2));

%   reconstructing the data set for one x at a time
for x=1:z(1)
    
    %   first matrix is A_x, nvc*nc, coil compression matrix
    %   second matrix is X_x, nc*ny, data matrix
    new(x,:,:)=sq(A(x,:,:))*sq(X(x,:,:));
    
end

%   place y-dim back to its place
new=permute(new, [1 3 2]);

%   inverse the fourier transform along x
new=fftshift(ifft(fftshift(new)));

disp('Coil Compression Completed')
toc

return

