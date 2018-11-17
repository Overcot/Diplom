function [ u ] = projection(u)
   global s t umin umax;
   for time=t
        for class=s
            if u(class, time) > umax
                u(class, time) = umax;
            elseif u(class, time) < umin
                u(class, time) = umin;
            end
        end
   end
end
