function value = saturate(value, limit, max)
    % Satura el valor dentro del rango [lower_limit, upper_limit]
    for i=1:3
       if (value(i) > limit) 
           value(i) = max;
       elseif (value(i) < -1*limit) 
           value(i) = -1*max;
       else
           value(i) = value(i);
       end
    end
end