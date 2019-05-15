#imports
import math
import numpy as np
import matplotlib.pyplot as plt
from ISA import *

def flightSim():
    print("A spacecraft is re-entering the atmosphere at 100km, don't panic.")
    t_interval = 10 #milliseconds

    rHeight = 325000*0.3048
    F_grav = [0,0]
    F_drag = [0,0]
    fp_angle = math.radians(-1.5)
    V_i_abs = 22500*.3048 
    V_i_x = V_i_abs*math.cos(fp_angle)
    V_i_y = V_i_abs*math.sin(fp_angle)
    V_SC = [V_i_x, V_i_y]
    speed_SC = ((V_SC[0]**2)+(V_SC[1]**2))**(1/2)
    a_SC = [0,0]
    

    G = 6.6741*10**(-11)
    E_m = 5.972*10**(24)
    E_r = 6371*10**(3)
    SC_m = 2662.8*0.4536
    SC_A = 2.81
    SC_CD = 1.6

    coord = [0,(rHeight+E_r)] #to centre of Earth
    coord_R = (coord[0]**2+coord[1]**2)**(1/2)
    D_earthsurface = round((coord_R - E_r)/1000, 1)

    E_p = -(G*E_m*SC_m)/(E_r+coord_R)
    E_k = 0.5*SC_m*speed_SC**(2)
    E_total = E_p+E_k

    coord_log = [[coord[0]/1000], [coord[1]/1000], [coord_R]]
    energy_log1 = [[0],[0], [E_total/10**9]]
    energy_log2 = [[],[]]

    #opening log file.
    positionLog = open("log_position.txt", "w+")

    timeMS=0
    while coord_R >= E_r:
        cHeight, cBoundaryName, cTemperature, cPressure, cDensity = heightAnalysis(D_earthsurface*1000)
        
        #radius and coordinate unit vector for gravity calculations
        coord_unit_vector = [coord[0]/(coord_R), coord[1]/(coord_R)]
        velocity_unit_vector = [V_SC[0]/speed_SC, V_SC[1]/speed_SC]

        #gravity calculations
        g = (G*E_m)/(E_r+coord_R)**(2)
        F_grav_magnitude = g*SC_m
        F_grav = [-F_grav_magnitude*coord_unit_vector[0], -F_grav_magnitude*coord_unit_vector[1]]

        #drag calculations
        F_drag_magnitude = (SC_CD*0.5*SC_A*cDensity*speed_SC**2)
        F_drag = [-1*F_drag_magnitude*velocity_unit_vector[0], -1*F_drag_magnitude*velocity_unit_vector[1]]
        if(cDensity > 0):
            tempTerminalVelocity = (g*SC_m/SC_CD*2/cDensity/SC_A)**(1/2)
            if(speed_SC <= abs(tempTerminalVelocity)):
                F_drag = [-F_grav[0], -F_grav[1]]

        #net_force
        F_net = [F_grav[0]+F_drag[0], F_grav[1]+F_drag[1]]

        #acceleration
        a_SC = [F_net[0]/SC_m, F_net[1]/SC_m]

        #coord update
        coord[0] += t_interval/1000*V_SC[0]
        coord[1] += t_interval/1000*V_SC[1]
        coord_R = (coord[0]**2+coord[1]**2)**(1/2)

        #velocity
        V_SC[0] += t_interval/1000*a_SC[0]
        V_SC[1] += t_interval/1000*a_SC[1]
        speed_SC = ((V_SC[0]**2)+(V_SC[1]**2))**(1/2)

        

        D_earthsurface = round((coord_R - E_r)/1000, 1)
        coord_log[0].append(coord[0]/1000)
        coord_log[1].append(coord[1]/1000)
        coord_log[2].append(coord_R/1000)

        E_p = -(G*E_m*SC_m)/(E_r+coord_R)
        E_k = 0.5*SC_m*speed_SC**(2)
        E_total = E_p+E_k
        E_change = E_total - energy_log1[2][-1]

        energy_log1[0].append(timeMS/1000/60)
        energy_log1[1].append(E_change/10**9)
        energy_log1[2].append(E_total/10**9)

        energy_log2[0].append(D_earthsurface)
        energy_log2[1].append(speed_SC)
        
        intervalPrint = 60
        if(timeMS % (1000*intervalPrint) == 0):

            #print("X-axis, coord[km]: " + str(round(coord[0]/1000,1)) + ", velocity[m/s]: " + str(round(V_SC[0],1)) + ", acceleration[m/s^2]: " + str(round(a_SC[0],5))   + ". F_grav [N]: " + str(round(F_grav[0],1))  + ". F_drag [N]: " + str(round(F_drag[0],1)) + "."   )
            #print("Y-axis, coord[km]: " + str(round(coord[1]/1000,1)) + ", velocity[m/s]: " + str(round(V_SC[1],3)) + ", acceleration[m/s^2]: " + str(round(a_SC[1],3))   + ". F_grav [N]: " + str(round(F_grav[1],0))  + ". F_drag [N]: " + str(round(F_drag[1],0)) + ". Net force[N]:"   + str(round(F_net[1],0)))
            print("Time[min]: " + str(round(timeMS/1000/60,1))+ ", DE[km]: " + str(round(D_earthsurface,1)) + ". Speed[m/s]:  " + str(round(speed_SC,1)) + ". X-Coord[km]: " + str(round(coord[0]/1000,1)) + ". Y-Coord[km]:" + str(round(coord[1]/1000,1)) + ". Total energy[GJ]: " + str(round(E_total/10**(9), 2)) + ".")
            
            
            logWrite = "\nTime[s]: " + str(round(timeMS/1000,1)) + ". DE[km]: " + str(round(D_earthsurface,1)) + ". Speed[m/s]:  " + str(round(speed_SC,1)) + ". X-Coord[km]: " + str(round(coord[0]/1000,1)) + ". Y-Coord[km]:" + str(round(coord[1]/1000,1))
            positionLog.write(logWrite)

        hourLimit = 30
        timeMSLimit = hourLimit*3600*1000
        if(timeMS > timeMSLimit):
            break
        

        timeMS += t_interval
    

    #plotting trajectory
    fig, ax = plt.subplots(2,2)

    root = fig.canvas._tkcanvas.winfo_toplevel()

    theta = np.linspace(0,2*np.pi,100)
    r = E_r/1000
    ax[0][0].plot((r*np.cos(theta)),(r*np.sin(theta)))
    ax[0][0].plot(coord_log[0], coord_log[1])
    ax[0][0].set_xlabel("X-Coord [km]")
    ax[0][0].set_ylabel("Y-Coord [km]")
    ax[0][0].axis("equal")

    ax[1][0].plot(energy_log1[0], energy_log1[1])
    ax[1][0].set_xlabel("Time [min]")
    ax[1][0].set_ylabel("Energy Change [GJ]")

    ax[1][1].plot(energy_log2[0], energy_log2[1])
    ax[1][1].set_xlabel("Distance to Earth's Surface [km]")
    ax[1][1].set_ylabel("Speed [m/s]")
    ax[1][1].invert_xaxis()

    plt.tight_layout()
    plt.show()
    root.mainloop()
        

#main loop
if __name__ == "__main__":
    
    generateBaseLayers()
    
    startUp = True

    #heightChooser analysis
    ###heightChooserAnalysis()

    #flightSim
    flightSim()

    



