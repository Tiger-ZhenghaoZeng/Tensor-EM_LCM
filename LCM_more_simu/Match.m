function [L1,w1,Theta1,v1] = Match(L,w,Theta,v)
% Match estimates and true parameters
n = length(L);
for i = 1:(n-1)
    k = i;
    min = norm(v(:,k)-Theta(:,i));
    for j = (i+1):n
        current = norm(v(:,j)-Theta(:,i));
        if (current < min)
            k = j;
            min = current;
        end
    end
    a = v(:,i);
    v(:,i) = v(:,k);
    v(:,k) = a;
    b = w(i);
    w(i) = w(k);
    w(k) = b;
end
L1 = L;
w1 = w;
Theta1 = Theta;
v1 = v;