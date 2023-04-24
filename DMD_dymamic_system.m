function [v_real_sig, Tu] = DMD_dynamic_system(A,B,v_real, v_ref,ts,Tu)

vp = A*v_real+B*v_ref - Tu;
v_real_sig = v_real + vp*ts;

end