# plot density of SS-pLBA model for a parameter setting
#
plot.density <- function(pars, RT.lims=c(0,4)){
  pls=make.parlists(pars)
  t=seq(RT.lims[1],RT.lims[2], by=0.01)
  par(mfrow=c(2,2))
  
  plot(t,dfun.1(t,pls$D1S1), type='l',
       xlim=c(-RT.lims[2],RT.lims[2]), 
       ylab='pdf', xlab='t (s)', main="normal GO left")
  lines(-t,dfun.2(t,pls$D1S1), type='l',xlim=c(-RT.lims[2],RT.lims[2]))

  plot(t,dfun.1(t,pls$D2S1), type='l',
       xlim=c(-RT.lims[2],RT.lims[2]), 
       ylab='pdf', xlab='t (s)', main="deprived GO left")
  lines(-t,dfun.2(t,pls$D2S1), type='l',xlim=c(-RT.lims[2],RT.lims[2]))
  
  plot(t,dfun.2(t,pls$D1S2), type='l',
       xlim=c(-RT.lims[2],RT.lims[2]), 
       ylab='pdf', xlab='t (s)', main="normal GO right")
  lines(-t,dfun.1(t,pls$D1S2), type='l',xlim=c(-RT.lims[2],RT.lims[2]))

  plot(t,dfun.2(t,pls$D2S2), type='l',
       xlim=c(-RT.lims[2],RT.lims[2]), 
       ylab='pdf', xlab='t (s)', main="deprived GO right")
  lines(-t,dfun.1(t,pls$D2S2), type='l',xlim=c(-RT.lims[2],RT.lims[2]))
}

#plot.model(pars)
#plot.model(untrans(fit$par))
