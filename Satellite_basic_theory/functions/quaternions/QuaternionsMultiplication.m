function qotimesq1 = QuaternionsMultiplication(q,q1)
%Function to calculate the Crassidis quaternion multiplication
%According to: 
%[1] Markley, F. Landis, and John L. Crassidis. Fundamentals of spacecraft 
%    attitude determination and control. Vol. 1286. New York, NY: Springer 
%    New York, 2014.

% Eq 2.84
psi = [-q(2) -q(3) -q(4);
        q(1)  q(4) -q(3);
       -q(4)  q(1)  q(2);
        q(3) -q(2)  q(1)];
      
qotimes = [q psi];
qotimesq1 = qotimes*q1; 
end

