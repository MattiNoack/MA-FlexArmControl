% INPUT: Torque 1
function u = torque1(t)

    if(t<10)
        u = 0.2*t;
    else
        u = 2;
    end

end