function Tdis = simplifiedDisturbances(t)
    f = 3.53E-4;
    a = 1.37E-7;
    phase = rand(3,1);
    off = 1.02E-7;
    Tdis = a*sin(2*pi*f*t+phase)+off;
end