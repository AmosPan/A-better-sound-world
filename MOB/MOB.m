%% Introduction
% It's the program of Multipole Orthogonal Beamforming by Xingjian Pan.
%% Contact
% Feel free to contact me: xj_pan@sjtu.edu.cn
%% 
% You need to input the cross-spectral matrix of array signals, positions
% of the microphone array (a Nm x 3 matrix), positions of the scanning
% plane (a Ng x 3 matrix) and the analysis frequency.
% the outputs include the beamforming outputs at each grid point and the
% positions of sources.
function [B,posi]=OB_multipole(CSM,mics,grids,f)
c0=343; % sound speed
k=2*pi*f/c0; % wavenumber
Ng=size(grids,1);Nm=size(mics,1);
% the number of grid points and the number of microphones
B=zeros(Ng,1); % initialization 
cen=mean(mics,1); % the array center
[U,S,V]=svd(CSM); % the singular value decomposition
% to compute the number of sources
eigen=diag(S); % the eigenvalues
m=eigen(1:end-1)./eigen(2:end); % the gradient of the eigenvalues
n=find(m==max(m)); % the number of sources
posi=zeros(n,3);
% to compute the steering vectors
h=zeros(Nm,Ng);
for i=1:Ng
    R_grid=ones(Nm,1)*grids(i,:)-mics;  
    r_grid=sqrt(R_grid(:,1).^2+R_grid(:,2).^2+R_grid(:,3).^2); 
    r_cen=norm(cen-grids(i,:));
    h(:,i)=1/Nm*r_grid/r_cen.*exp(-2*1i*k*(r_grid-r_cen));
end
% The beamforming procedure
for j=1:n
    u=U(:,j);
    factor=sum(abs(sqrt(S(j,j))*u))/Nm;% the correction for the source strength
    u=u./abs(u); % the normalization procedure
    u=u.^2; % the square procedure
    C_sub=u*u'; % the modified cross-spectral matrix
    B_sub=zeros(Ng,1);
    for i=1:Ng
      B_sub(i)=factor^2*real(h(:,i)'*C_sub*h(:,i)); % The beamforming procedure.
    end
    num=find(B_sub==max(B_sub));num=num(1);
    grid_max=grids(num,:); 
    posi(j,:)=grid_max;
    R=grids-ones(Ng,1)*grid_max;
    r=sqrt(R(:,1).^2+R(:,2).^2+R(:,3).^2);
    Q=max(B_sub)*10.^(-500*r.^2); % the clean procedure for a higher resolution
    B=B+Q; % The final source map is a linear combination of beamforming outputs of different sources.
end