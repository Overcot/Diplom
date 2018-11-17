function f=AnnaModel(ab)
global recruitment ssb maxSSB
f = 0;
for t = 1:length(recruitment)
    if ssb(t)/maxSSB > ab(3)
        f = f + (log(recruitment(t)) - log(ssb(t)) - ab(1) + ab(2)*ssb(t) ).^2;
    else
        f = f + (log(recruitment(t)) - log(ssb(t)) - ab(1) + ab(2)*ssb(t) - (ssb(t)/maxSSB-ab(3)) ).^2;
    end
end
