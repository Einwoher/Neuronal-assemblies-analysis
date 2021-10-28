function z = baseline(y, lambda, p) %from Eilers&Boelens, 2005)
% Estimate baseline with asymmetric least squares 
m = length(y); 
D = diff(speye(m), 2); 
w = ones(m, 1); 
for it = 1:100
W = spdiags(w, 0, m, m);
D0=D';
C = chol(W+lambda*D0*D); 
C0=C';
z = C\(C0\(w .* y)); 
w = p * (y > z) + (1 - p) * (y < z);
end