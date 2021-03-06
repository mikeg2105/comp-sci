function [psisum]=beat_wave(t,amp1,amp2, om, fac, ph)
% creating two superimposed cos functions...

   psisum=psidef(t, amp1, om, ph)+psidef(t, amp2, fac*om, ph);
end

function [psi]=psidef(t, A, omega_0, phi)
    psi=A*cos(2*pi*omega_0*t+2*pi*2*pi*phi/360);
end

function [psi]=psidefv(t, A, k, r, omega_0, phi)
    psi=A*cos(k*r+2*pi*omega_0*t+2*pi*2*pi*phi/360);
end