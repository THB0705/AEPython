#imports
import math
import numpy as np
import matplotlib.pyplot as plt


#constants
g_0 = 9.80665
R = 287
T_0 = 288.15
p_0 = 101325
a_layer = [-6.5, 0, 1, 2.8, 0, -2.8, -2, 0, 0]
a_boundaries = [11, 20, 32, 47, 51, 71, 86, 100, 1000000000]
a_boundaries_names = ["Troposphere", "Tropopause", "Stratosphere", "Stratosphere", "Stratopause", "Mesosphere", "Mesosphere", "preThermoSphere", "space"]
T_boundaries = [T_0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
P_boundaries = [p_0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#booleans


#Generating Baselayers
def generateBaseLayers():
    layerG = 0
    while layerG < 8:
        if(layerG == 0):
            T_boundaries[layerG+1] = T_boundaries[layerG] + a_layer[layerG]*a_boundaries[layerG]
        elif(layerG != 0):
            T_boundaries[layerG+1] = T_boundaries[layerG] + a_layer[layerG]*(a_boundaries[layerG]-a_boundaries[layerG-1])

        if(a_layer[layerG] == 0):
            P_boundaries[layerG+1] = isothermalLayer(P_boundaries[layerG], T_boundaries[layerG], 1000*(a_boundaries[layerG]-a_boundaries[layerG-1]))
        elif(a_layer[layerG] != 0):
            P_boundaries[layerG+1] = normalLayer(P_boundaries[layerG], T_boundaries[layerG], T_boundaries[layerG+1], a_layer[layerG]/1000)
        #print("Layer No.: " + str(layerG) + ". T_0: " + str(round(T_boundaries[layerG],2)) + " K. P_0: " + str(round(P_boundaries[layerG],2)) + " Pa.")

        layerG += 1

#input loop
def menuInput(startUp):
    if startUp == True:
        print("    **** ISA Calculator ****")
    
    height = 0
    unit = "meter"

    while True:
        unitIn = input("In what units will you enter your height?\n1. Meters\n2. Feet\n3. Flight Levels (FL)\n")
        try:
            unitInTemp = int(unitIn)
            if(unitInTemp > 0 and unitInTemp < 4):
                if(unitInTemp == 1):
                    unit = "meter"
                    break
                elif(unitInTemp == 2):
                    unit = "feet"
                    break
                elif(unitInTemp == 3):
                    unit = "FL"
                    break
            else:
                print("You did not enter a number between 1 and 3.")
        except ValueError:
            print("You did not enter a number. Please try again.")

    while True:
        heightInput = input("\nEnter height: ")
        try:
            heightConversion = float(heightInput)
            if(unit == "meter"):
                height = heightConversion
            elif(unit == "feet"):
                height = heightConversion*0.3048
            elif(unit == "FL"):
                height = heightConversion*0.3048*100
            if(height < 0 or height > 100000):
                print("The range of heights we can calculate atmospheres for is 0m to 100km. Enter a number that fits the range.")
            else:
                break
        except ValueError:
            print("Error, make sure you enter a number. Try again please.")
    
    return height

#specific Data Generation

def heightAnalysis(height):
    x = 0

    while x in range(len(a_boundaries)):
        boundary = a_boundaries[x]*1000
        if(height <= boundary):
            layer = x
            break
        x += 1   

    pressureH = 0
    rhoH = 0
    temperatureH = 1
    if(height<=100*10**3):

        if(layer == 0):
            temperatureH = T_boundaries[layer] + a_layer[layer]/1000*height
            if(a_layer[layer] == 0):
                pressureH = isothermalLayer(P_boundaries[layer], T_boundaries[layer], 1000*a_boundaries[layer])
            elif(a_layer[layer] != 0):
                pressureH = normalLayer(P_boundaries[layer], T_boundaries[layer], temperatureH, a_layer[layer]/1000)
        elif(layer != 0):
            temperatureH = T_boundaries[layer] + a_layer[layer]/1000*(height-a_boundaries[layer-1]*1000)
            if(a_layer[layer] == 0):
                pressureH = isothermalLayer(P_boundaries[layer], T_boundaries[layer], (height-1000*a_boundaries[layer-1]))
            elif(a_layer[layer] != 0):
                pressureH = normalLayer(P_boundaries[layer], T_boundaries[layer], temperatureH, a_layer[layer]/1000)
    
        rhoH = pressureH/R/temperatureH

    return height, a_boundaries_names[layer], temperatureH, pressureH, rhoH

#specific functions
def isothermalLayer(p_0, T, delta_H):
    p_H = p_0 * math.exp(-g_0/R/T*delta_H)
    return p_H

def normalLayer(p_0, T_0, T_H, a):
    p_H = p_0 * (T_H/T_0)**(-g_0/(a)/R)
    return p_H

def printStatementGenerator(height, boundaryName, temperature, pressure, density):
    printStatement = "\n***\n" + "Height:      " + str(height) + " m.\n" + "Layer:       " + boundaryName + ".\n"+ "Temperature: " + str(round(temperature)) + " K.\n"+ "Pressure:    " + str(round(pressure,2)) + " Pa.\n"+ "Density:     " + str(round(density,4)) + " kg/m^3."
    return printStatement    

def heightChooserAnalysis():
    height = menuInput(startUp = True)

    cHeight, cBoundaryName, cTemperature, cPressure, cDensity = heightAnalysis(height)

    finalPrint = printStatementGenerator(cHeight, cBoundaryName, cTemperature, cPressure, cDensity)

    print(finalPrint)

def flightSim():
    print("A spacecraft is re-entering the atmosphere at 100km, don't panic.")
    t_interval = 10 #milliseconds

    rHeight = 325000*0.3048
    F_grav = [0,0]
    F_drag = [0,0]
    fp_angle = math.radians(-10)
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

    coord_log = [[coord[0]/1000], [coord[1]/1000], [coord_R]]

    #flightSimLoop

    timeMS=0
    while coord_R >= E_r:
        cHeight, cBoundaryName, cTemperature, cPressure, cDensity = heightAnalysis(D_earthsurface*1000)
        
        #radius and coordinate unit vector for gravity calculations
        coord_unit_vector = [coord[0]/(coord_R), coord[1]/(coord_R)]
        velocity_unit_vector = [V_SC[0]/speed_SC, V_SC[1]/speed_SC]

        #gravity calculations
        F_grav_magnitude = (G*E_m*SC_m)/(E_r+coord_R)**(2)
        F_grav = [-F_grav_magnitude*coord_unit_vector[0], -F_grav_magnitude*coord_unit_vector[1]]

        #drag calculations
        F_drag_magnitude = (SC_CD*0.5*SC_A*cDensity*speed_SC**2)
        F_drag = [-1*F_drag_magnitude*velocity_unit_vector[0], -1*F_drag_magnitude*velocity_unit_vector[1]]
        #if(abs(F_drag[0]) > abs(F_grav[0])):
            #F_drag = 

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

        
        
        if(timeMS % (1000*60) == 0):

            #print("X-axis, coord[km]: " + str(round(coord[0]/1000,1)) + ", velocity[m/s]: " + str(round(V_SC[0],1)) + ", acceleration[m/s^2]: " + str(round(a_SC[0],5))   + ". F_grav [N]: " + str(round(F_grav[0],1))  + ". F_drag [N]: " + str(round(F_drag[0],1)) + "."   )
            #print("Y-axis, coord[km]: " + str(round(coord[1]/1000,1)) + ", velocity[m/s]: " + str(round(V_SC[1],3)) + ", acceleration[m/s^2]: " + str(round(a_SC[1],3))   + ". F_grav [N]: " + str(round(F_grav[1],0))  + ". F_drag [N]: " + str(round(F_drag[1],0)) + ".   :"   + str(cDensity))
            print("Time[min]: " + str(round(timeMS/1000/60,1))+ ", DE[km]: " + str(round(D_earthsurface,1)) + ". Speed[m/s]:  " + str(round(speed_SC,1)) + ". X-Coord[km]: " + str(round(coord[0]/1000,1)) + ". Y-Coord[km]:" + str(round(coord[1]/1000,1)) + ".")

        timeMS += t_interval
    

    #plotting trajectory
    fig, ax = plt.subplots()

    root = fig.canvas._tkcanvas.winfo_toplevel()

    theta = np.linspace(0,2*np.pi,100)
    r = E_r/1000
    ax.plot((r*np.cos(theta)),(r*np.sin(theta)))

    ax.plot(coord_log[0], coord_log[1])
    plt.xlabel("X-Coord [km]")
    plt.ylabel("Y-Coord [km]")
    
    ax.axis("equal")

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

    



