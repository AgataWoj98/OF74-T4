import numpy as np
import matplotlib.pyplot as plt 
import time
import sys


def acceleration(v):
    return np.array([0,-g])-f*np.linalg.norm(v)*v


d = 6 #odległosc tenisistki od siatki w momencie uderzenia
hs = 1 #wysokosc siatki
a = 12 #odległosc rogu kortu od siatki 
b = 11 #szerokosc kortu
m = 0.058 #masa piłki
rad = 0.0325 #pomień piłki
ro = 1.2 #gęstosć powietrza
c = 0.45 #współczynnik oporu powietrza
g = 9.81 #przyspieszenie ziemskie 

f = c*np.pi*rad**2*ro/2/m #siła oporu F=f*m*v**2
l0 = np.sqrt((a+d)**2+(b/2)**2) #długosc rzutu 
ls = l0*d/(a+d) #odległosc tenisitki od siatki

alphas = np.linspace(0.220, 0.223, 500) #skanowane kąty prędkosci początkowej wzgl poziomu
v0s = np.linspace(22.8, 23.1, 500) #skanowane wartosci prędkosci początkowej

h0 = 0 #początkowa wysokosc piłki
print(alphas, v0s)


dt = 0.0001
T_max = 20
T = np.arange(0, T_max, dt) #czas propagacji
T_throw = np.zeros((len(v0s), len(alphas)))

l_err = 0.02 #margines błędu zasięgu rzutu


#wykres toru ruchu z siatką (czarne) i kortem (zielone)
plt.figure()
plt.xlabel("x")
plt.ylabel("y")
plt.vlines(ls, 0, hs, colors="black")
plt.hlines(0, 0, l0, colors="green")


t_min = T_max
for id_v0, v0 in enumerate(v0s):
    print(id_v0)
    start = time.time()
    for id_alpha, alpha in enumerate(alphas):
        
        #inicjacja list położeń i prędkosci i ich wartosci początkowych 
        r = []
        v = []
        r.append(np.array([0,h0]))
        v.append(np.array([v0*np.cos(alpha), v0*np.sin(alpha)])) 
        
        over = 0 #czy piłka przeszła przez siatkę
        
        #propagacja algorytmem Eulera
        for i in range(0, len(T)-1):
            t = T[i]
            acc = acceleration(v[i])
            r.append(r[i] + v[i]*dt + acc*dt**2/2)
            v.append(v[i] + acc*dt)
            

            
            if(r[i+1][0]-l0>l_err): #czy piłka wyrzucona na out?
                break
                
            if(r[i+1][1]>hs and r[i+1][0]>=ls): #czy piłka przeszła przez siatkę?
                over = 1
                
            if(over==0 and r[i+1][0]>=ls): #piłka nie przeszła przez siatkę
                break
            
            if(r[i+1][1]<=0): #piłka dotknęła ziemi za siatką, po minimalnym czasie na zasięgu l0 +- l_err
                if(over==1 and np.abs(r[i+1][0]-l0)<=l_err and t<t_min): 
                    T_throw[id_v0][id_alpha] = t
                    plt.scatter(np.array(r)[:,0], np.array(r)[:,1], s=1)
                    t_min = t
                    v0_opt = v0
                    alpha_opt = alpha
                break
        #np.save("trajectory"+str(v0)+"_"+str(alpha)+".npy", r) dobrze jest zapisywać dane na wszelki wypadek
    print("time for one v0 ", time.time()-start)

        

print(T_throw)
print("Minimalny czas rzutu=", t_min, "s dla v0=$ ", v0_opt, "m/s i alpha=", alpha_opt)

plt.figure()
plt.imshow(T_throw!=0)

plt.figure()
nticks = 6
ticks = [int(i) for i in range(len(alphas)) if(i%(len(alphas)//nticks)==0)]
plt.title(r"$t_{min}$="+str(np.round(t_min,3))+r"s $v_0$="+str(np.round(v0_opt,3))+r"$\frac{m}{s}$ $\alpha$="+str(np.round(alpha_opt,3))+" dt="+str(dt)+r"s $l_{err}$="+str(l_err)+"m")
plt.imshow(T_throw, vmin=0.975, vmax=0.98)
plt.xticks(ticks, labels=np.round(alphas[ticks]*180/np.pi,1))
plt.yticks(ticks, labels=np.round(v0s[ticks],1))
plt.xlabel(r"$\alpha$ [$^o$]")
plt.ylabel(r"$v_0$ [$\frac{m}{s}$]")
plt.colorbar()

#plt.savefig("final_500.pdf")



#porównanie toru z parabolą (dla f=0 scisłe)
"""
plt.figure()
plt.xlabel("x")
plt.ylabel("(Euler-parabola)/parabola [%]")
plt.plot(np.array(r)[:,0], (np.array(r)[:,1] + g/2/v0**2 * r[:,0]**2)/(g/2/v0**2 * r[:,0]**2)*100)
"""





