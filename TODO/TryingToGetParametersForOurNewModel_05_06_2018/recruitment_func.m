function f = recruitment_func(ab)
global recruitment ssb
f = 0;
for t = 1:length(recruitment)
    f = f + (log(recruitment(t)) - log(ab(1)*ssb(t)) - (-ab(2).*ssb(t)) ).^2;
end
