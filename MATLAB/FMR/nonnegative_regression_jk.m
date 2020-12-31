function [ beta,exitflag ] = nonnegative_regression_jk( X,y,varargin )
%Weighted linear regression of y on X with weights W
%   Detailed explanation goes here

p=inputParser;

addRequired(p,'X',@(x)ismatrix(x) & isreal(x))
[nn,pp]=size(X);
if pp>nn
    warning('This method of solving nonnegative least squares is inefficient for underdetermined systems')
end
addRequired(p,'y',@(x)isvector(x) && size(x,1)==nn)
addOptional(p,'W',ones(nn,1),@(x)isvector(x) && size(x,1)==nn)%
addOptional(p,'JackknifeBlocks',[],@(x)iscell(x))
addOptional(p,'NoJackknifeBlocks',100,@(x)isscalar(x))
parse(p,X,y,varargin{:} )

JackknifeBlocks=p.Results.JackknifeBlocks;
NoJackknifeBlocks=p.Results.NoJackknifeBlocks;
W=p.Results.W;

if isempty(JackknifeBlocks)
    JackknifeBlocks=cell(NoJackknifeBlocks,1);
    blocksize=ceil(nn/NoJackknifeBlocks);
    for ii=1:NoJackknifeBlocks
        JackknifeBlocks{ii}=(ii-1)*blocksize+1:min(nn,ii*blocksize);
    end
else
    NoJackknifeBlocks=length(JackknifeBlocks);
end

alpha=zeros(pp,NoJackknifeBlocks);
Sigma=zeros(pp,pp,NoJackknifeBlocks);
beta=zeros(NoJackknifeBlocks,pp);
exitflag=zeros(NoJackknifeBlocks,1);
for ii=1:NoJackknifeBlocks
    alpha(:,ii)=X(JackknifeBlocks{ii},:)'*diag(sparse(W(JackknifeBlocks{ii})))*y(JackknifeBlocks{ii});
    Sigma(:,:,ii)=X(JackknifeBlocks{ii},:)'*diag(sparse(W(JackknifeBlocks{ii})))*X(JackknifeBlocks{ii},:);
end
options=optimoptions('lsqlin','Display','None');
for ii=1:NoJackknifeBlocks
    incl=[1:ii-1,ii+1:NoJackknifeBlocks];
    [beta(ii,:),~,~,exitflag(ii)]=lsqlin(sum(Sigma(:,:,incl),3)/nn,sum(alpha(:,incl),2)/nn,-eye(pp),zeros(pp,1),[],[],[],[],[],options);
end
%beta=(X'*diag(sparse(W))*X)\(X'*diag(sparse(W))*y);


end

