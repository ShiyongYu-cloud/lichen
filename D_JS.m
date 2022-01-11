function dist = D_JS(P,Q)
%% function for calculating the Jensen-Shannon divergence of two PDFs
% P and Q are automatically normalized to have the sum of one
%%
P = P(:);
Q = Q(:);
% normalizing the P and Q
P = P./sum(P);
Q = Q./sum(Q);
M = (P + Q)/2;
d1 = D_KL(P,M);
d2 = D_KL(Q,M);
dist = (d1 + d2)/2;
return;
%%
function dist = D_KL(P1,P2)
%% function for calculating the Kullback-Leibler divergence of two PDFs
X = log2(P1./P2);
X(isnan(X)) = 0; % resolving the case when P(i)==0 & Q(i)==0
Y = P1.*X;
Y(isnan(Y)) = 0; % resolving the case when P(i)==0   
dist = sum(Y);
return;