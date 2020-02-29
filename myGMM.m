function y = myGMM(ran,mu,sigma)
y = zeros(1,length(ran));
for i = 1:length(mu)
    y = y + normpdf(ran,mu(i),sigma);
end
y = y./sum(y);