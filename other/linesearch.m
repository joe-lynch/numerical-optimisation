function alpha = linesearch(x,fx,gx,s,theta,f)

% Backtracking line search algorithm

% Initialise alpha and set function increment for sufficient descent 
% condition

alpha = 1;
df = theta*gx'*s;

% Half alpha until the sufficient descent condition is satisfied; stop when
% alpha = 2^{-30} < 10^{-9}

for i=1:30

    if f(x+alpha*s) - fx < alpha*df 
        break;
    end
    
    alpha = alpha/2;
        
end

end

