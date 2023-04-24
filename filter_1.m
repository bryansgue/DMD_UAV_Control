function [learning, F_filter, F_memory] = filter_1(F, F_filter, F_memory, A, B)

% Learning algorithm
l1_k = -B*F_filter + A*F_memory;


% Vector of learning parameters
learning = [l1_k];

% Update values memories learning 1
F_filter = update_l_memories(learning, F_filter);
F_memory = update_e_memories(F, F_memory);

end

function l1_memories = update_l_memories(l, l1_memories)
for k = length(l1_memories):-1:2
    
    l1_memories(k) = l1_memories(k-1);
    
end
l1_memories(1) = l(1);

end

function e1_memories = update_e_memories(xe, e1_memories)
for k = length(e1_memories):-1:2
    
    e1_memories(k) = e1_memories(k-1);
   
end
e1_memories(1) = xe(1);

end