function [ value_f_bs ] = Forward_BS( strike, tenor, spot, rateclp, rateusd, vol )

[ value_call ] = ValueBS( strike, tenor, spot, rateclp, rateusd, vol, 1);
[ value_put ] = ValueBS( strike, tenor, spot, rateclp, rateusd, vol, -1);

value_f_bs=value_call-value_put;

end

