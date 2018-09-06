function f = ricker(ab, xdata)

f = log(ab(1).*xdata.*exp(-ab(2).*xdata));
end

