function x = timestep_forwardEuler(ueval,veval,weval,X,dt)

x = X + dt*[ueval,veval,weval];

end