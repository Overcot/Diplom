function f = ricker(ab)
global ssb recruitment
%f = log(ab(1).*xdata.*exp(-ab(2).*xdata));
f = log(ab(1).*ssb.*exp(-ab(2).*ssb)) - log(recruitment);
end

