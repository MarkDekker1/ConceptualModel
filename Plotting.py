#%%
tvec=np.linspace(0,tmax,np.int(tmax/Dt))

plt.figure(num=None, figsize=(8,4),dpi=150, facecolor='w', edgecolor='k')
plt.plot(tvec,m1.T_eq_ocean, '-',linewidth=2)
plt.plot(tvec,m1.T_pole_ocean, '-',linewidth=2)
plt.plot(tvec,m1.T_eq_atm, '-',linewidth=2)
plt.plot(tvec,m1.T_pole_atm, '-',linewidth=2)
plt.ylabel(r'Quantity',fontsize=15)
plt.xlabel('Time [s]',fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.ylim([200,320])
plt.legend(['Eq_ocean','Pole_ocean','Eq_atm','Pole_atm'])

#%%Growth rates
G_pole_ocean=[0]
G_eq_ocean=[0]
G_pole_atm=[0]
G_eq_atm=[0]
for i in range(1,np.int(tmax/Dt)):
    
    G_pole_ocean.append(m1.T_pole_ocean[i]-m1.T_pole_ocean[i-1])
    G_eq_ocean.append(m1.T_eq_ocean[i]-m1.T_eq_ocean[i-1])
    G_pole_atm.append(m1.T_pole_atm[i]-m1.T_pole_atm[i-1])
    G_eq_atm.append(m1.T_eq_atm[i]-m1.T_eq_atm[i-1])


plt.figure(num=None, figsize=(8,4),dpi=150, facecolor='w', edgecolor='k')
plt.semilogy(tvec,np.abs(G_eq_ocean), '-',linewidth=2)
plt.semilogy(tvec,np.abs(G_pole_ocean), '-',linewidth=2)
plt.semilogy(tvec,np.abs(G_eq_atm), '-',linewidth=2)
plt.semilogy(tvec,np.abs(G_pole_atm), '-',linewidth=2)
plt.ylabel(r'Quantity',fontsize=15)
plt.xlabel('Time [s]',fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.legend(['Eq_ocean','Pole_ocean','Eq_atm','Pole_atm'])

#%%Chemicals
plt.figure(num=None, figsize=(8,4),dpi=150, facecolor='w', edgecolor='k')
plt.semilogy(tvec,m1.Ca, '-',linewidth=2)
plt.semilogy(tvec,m1.Co, '-',linewidth=2)
plt.semilogy(tvec,m1.A, '-',linewidth=2)
plt.semilogy(tvec,m1.P, '-',linewidth=2)
plt.semilogy(tvec,m1.CO3, '-',linewidth=2)
plt.ylabel(r'Quantity',fontsize=15)
plt.xlabel('Time [s]',fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.legend(['Ca','Co','A','P','CO3'])
#%% Animation
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, np.int(tmax/Dt)), ylim=(250,350))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = tvec[0:i]
    y = m1.T_eq_ocean[0:i]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=250, interval=10, blit=True)

#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()