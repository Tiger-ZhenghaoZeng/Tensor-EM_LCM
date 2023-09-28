function [p1,w1,Theta1,v1] = Match(p,w,Theta,v)
% Match estimates and true parameters (required since latent classes are identified up to permutations)
% p, Theta: The proportion vector and item parameters as reference
% w, v: The proportion vector and item parameters to match with p, Theta
n = length(p);
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
p1 = p;
w1 = w;
Theta1 = Theta;
v1 = v;