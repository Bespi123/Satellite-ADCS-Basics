function [conj_q] = QuaternionsConjugate(q)
%Function to calculate the conjugate function
%According to: 
%[1] Markley, F. Landis, and John L. Crassidis. Fundamentals of spacecraft 
%    attitude determination and control. Vol. 1286. New York, NY: Springer 
%    New York, 2014.
conj_q = [q(1),-q(2),-q(3),-q(4)]';  %eq(2.91)
end

