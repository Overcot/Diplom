function f = recruitment_func(ab)
global recruitment ssb
f = 0;
for t = 1:53
    %f = f + (log(recruitment(t))-log(ssb(t)).*(ab(1)-ab(2).*ssb(t))).^2;
    f = f + (log(recruitment(t)) - (log(ab(1) .* ssb(t)) + log(-ab(2).*ssb(t)))).^2;
end
