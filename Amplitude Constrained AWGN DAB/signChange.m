
function [indices] = signChange (v)
t = v(1:(length(v)-1)).* v(2:length(v));
FirstSignPositive = v(1:(length(v)-1)) >= 0;
indices = [];
for i = 1:(length(t))
    if ((t(i) <0 && FirstSignPositive(i)) || (v(i) == 0 && v(i-1)>0))
        indices = [indices i];
    end
end

