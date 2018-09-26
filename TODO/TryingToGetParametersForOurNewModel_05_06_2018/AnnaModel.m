function f=AnnaModel(ab)
global recruitment ssb
f = 0;
for t = 1:length(recruitment)
    if ssb(t) > ab(3)
        f = f + (log(recruitment(t)) - log(ssb(t)*exp(ab(1) - ab(2)*ssb(t)))  ).^2;
    else
        f = f + (log(recruitment(t)) - log(ssb(t)*exp(ab(1) - ab(2)*ssb(t))*exp((ab(3)-ssb(t))) )  ).^2;
    end
end
