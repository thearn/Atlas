% This function compares the polynomial approximation of the complete ellipic
% integrals given in Abramowitz and Stegun with Matlab's implementation.

% The results are idential if m is set to m-1 in the polynomial
% approximation.

% testRing.m confirms that using the function as given by Matlab gives the
% correct result.

m = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];    

a0 = 1.3862944;
a1 = 0.1119723;
a2 = 0.0725296;
b0 = 0.5;
b1 = 0.1213478;
b2 = 0.0288729;
K = (a0+a1*m+a2*m.^2) + (b0+b1*m+b2*m.^2).*log(ones(size(m))./m);
[K2,E2] = ellipke(m);

subplot(1,2,1);
plot(m,K,m,K2);

a1 = 0.4630151;
a2 = 0.1077812;
b1 = 0.2452727;
b2 = 0.0412496;
E = (1+a1*m+a2*m.^2) + (b1*m+b2*m.^2).*log(ones(size(m))./m);

subplot(1,2,2);
plot(m,E,m,E2);