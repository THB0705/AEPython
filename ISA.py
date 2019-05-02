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
        heightInput = input("\nEnter height[m]: ")
        try:
            heightConversion = float(heightInput)
            if(unit == "meter"):
                height = heightConversion
            elif(unit == "feet"):
                height = heightConversion*0.3048
            elif(unit == "FL"):
                height = heightConversion*0.3048*100
            if(height < 0 or height > 86000):
                print("The range of heights we can calculate atmospheres for is 0m to 86000m. Enter a number that fits the range.")
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
    temperatureH = 0

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

    print(
        "\n***\n"+
        "Height:      " + str(height) + " m.\n" + 
        "Layer:       " + a_boundaries_names[layer] + ".\n"+
        "Temperature: " + str(round(temperatureH,2)) + " K.\n"+
        "Pressure:    " + str(round(pressureH,2)) + " Pa.\n"+
        "Density:     " + str(round(rhoH,4)) + " kg/m^3.")

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


