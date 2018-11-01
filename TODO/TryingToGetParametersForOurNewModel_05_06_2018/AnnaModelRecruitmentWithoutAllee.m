function recruitment=AnnaModelRecruitmentWithoutAllee(ab)
global ssb maxSSB lowerSSB higherSSB
recruitment = [];
for t = 1:length(ssb)
    value = log(ssb(t)) + ab(1) - ab(2)*ssb(t);
    recruitment = [recruitment; value];
end