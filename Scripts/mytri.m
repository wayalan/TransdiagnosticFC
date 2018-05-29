function x = mytri(A)
% 90*90*nsubj --->     nsubj * 4005
% put elements in upper trianglar matrix into a vector, 
[a,~,c] = size(A);
x = zeros(c, (a^2-a)/2);
t = 1;
for i = 1 : a-1
    for j = i+1 : a
        x(:,t) = squeeze(A(i,j,:));
        t = t + 1;
    end
end