function recruitment=AnnaModelRecruitment(ab)
global ssb maxSSB lowerSSB higherSSB
recruitment = [];
for t = 1:length(ssb)
    if ssb(t)/maxSSB > ab(3)
        value = log(ssb(t)) + ab(1) - ab(2)*ssb(t);
        lowerSSB = lowerSSB + 1;
    else
        value = log(ssb(t)) + ab(1) - ab(2)*ssb(t) + (ssb(t)/maxSSB-ab(3));
        higherSSB = higherSSB + 1;
    end
    recruitment = [recruitment; value];
end