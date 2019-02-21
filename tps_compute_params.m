function [affine non_affine]=tps_compute_params(Q1,Q2,R,Y1,K,lambda,varargin)
% Y is the tagent points n*2 or n*3
% param=[affine:non_affine]
% input a column vector as a weight vector for those points
Y=tps_n_c(Y1);
n=size(K,1);
if isempty(varargin)
    w=diag(ones(n,1));
else
    if length(varargin)==1
    w=diag(varargin{1});
    else
        error('only accept one more paramaters');
    end
end
if size(w)~=[n,n]
    error('weight vector need to be same length as the # of points');
    
end
M=K+n*lambda*inv(w^2);
non_affine=(Q2/(Q2'*M*Q2))*Q2'*Y;
affine=(R\Q1')*(Y-M*non_affine);
