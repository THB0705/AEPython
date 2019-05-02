#imports
import math


#constants
g_0 = 9.80665
R = 287
T_0 = 288.15
p_0 = 101325
a_layer = [-6.5, 0, 1, 2.8, 0, -2.8, -2]
a_boundaries = [11, 20, 32, 47, 51, 71, 86]
a_boundaries_names = ["Troposphere", "Tropopause", "Stratosphere", "Stratosphere", "Stratopause", "Mesosphere", "Mesosphere"]
T_boundaries = [T_0, 0, 0, 0, 0, 0, 0, 0]
P_boundaries = [p_0, 0, 0, 0, 0, 0, 0, 0]

#booleans


#Generating Baselayers
def generateBaseLayers():
    layerG = 0
    while layerG < 7:
        if(layerG == 0):
            T_boundaries[layerG+1] = T_boundaries[layerG] + a_layer[layerG]*a_boundaries[layerG]
        else:
            T_boundaries[layerG+1] = T_boundaries[layerG] + a_layer[layerG]*(a_boundaries[layerG]-a_boundaries[layerG-1])

        if(a_layer[layerG] == 0):
            P_boundaries[layerG+1] = isothermalLayer(P_boundaries[layerG], T_boundaries[layerG], 1000*(a_boundaries[layerG+1]-a_boundaries[layerG]))
        else:
            P_boundaries[layerG+1] = normalLayer(P_boundaries[layerG], T_boundaries[layerG], T_boundaries[layerG+1], a_layer[layerG]/1000)

        #print("Layer No.: " + str(layerG) + ". T_0: " + str(round(T_boundaries[layerG],2)) + " K. P_0: " + str(round(P_boundaries[layerG],2)) + " Pa.")

        layerG += 1

#input loop
def menuInput(startUp):
    if startUp == True:
        print("    **** ISA Calculator Troposphere ****")
    
    height = 0

    while True:
        heightInput = input("\nEnter height: ")
        try:
            height = float(heightInput)
            break
        except ValueError:
            print("Error, you did not enter a number. Try again.")
    
    return height

#specific Data Generation, combine following two functions  

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
    temperatureH = 0

    if(layer == 0):
        temperatureH = T_boundaries[layer] + a_layer[layer]/1000*height
        if(a_layer[layer] == 0):
            pressureH = isothermalLayer(P_boundaries[layer], T_boundaries[layer], 1000*a_boundaries[layer])
        else:
            pressureH = normalLayer(P_boundaries[layer], T_boundaries[layer], temperatureH, a_layer[layer]/1000)
    else:
        temperatureH = T_boundaries[layer] + a_layer[layer]/1000*(height-a_boundaries[layer-1]*1000)
        if(a_layer[layer] == 0):
            pressureH = isothermalLayer(P_boundaries[layer], T_boundaries[layer], (height-1000*a_boundaries[layer-1]))
        else:
            pressureH = normalLayer(P_boundaries[layer], T_boundaries[layer], temperatureH, a_layer[layer]/1000)
    
    rhoH = pressureH/R/temperatureH

    print("Layer: " + a_boundaries_names[layer] + ".\nHeight: " + str(height) + " m.\nTemperature: " + str(round(temperatureH,2)) + " K.\nPressure: " + str(round(pressureH,2)) + " Pa.\nDensity: " + str(round(rhoH,2)) + " kg/m^3.")

#specific functions
def isothermalLayer(p_0, T, delta_H):
    p_H = p_0 * math.exp(-g_0/R/T*delta_H)
    return p_H

def normalLayer(p_0, T_0, T_H, a):
    p_H = p_0 * (T_H/T_0)**(-g_0/(a)/R)
    return p_H


#main loop
if __name__ == "__main__":
    
    generateBaseLayers()
    
    startUp = True

    height = menuInput(startUp = True)

    heightAnalysis(height)


