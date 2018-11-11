function [ y ] = error_fxn( x )
%error_fxn - compute the error function of input x

err = @(t) exp(-(1000*t).^2);
coder.extrinsic('integral');
y = 2/sqrt(pi)*integral(err,0,x);

end

