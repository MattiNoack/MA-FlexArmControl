% INPUT: Torque 2
function u = torque2(t)

    if(t<10)
        u = 0.2*t;
    else
        u = 0;%-0.5*(t-10);
    end

end