#imports



#constants
g_0 = 9.80665
R = 287
T_0 = 288.15
p_0 = 101325
a_layer = [-6.5, 0, 1, 2.8, 0, -2.8, -2]
a_boundaries = [11, 20, 32, 47, 51, 71, 86]
a_boundaries_names = ["Troposphere", "Tropopause", "Stratosphere", "Stratosphere", "Stratopause", "Mesosphere", "Mesosphere"]

#booleans



def menuInput(startUp):
    if startUp == True:
        print("    **** ISA Calculator Troposphere ****")
    
    height = 0

    while True:
        heightInput = input("Enter height (meters): ")
        try:
            height = float(heightInput)
            break
        except ValueError:
            print("Error, make sure to enter your height as a number.")
    
    return height

def layerDetermination(height):

    x = 0

    while x in range(len(a_boundaries)):
        x += 1
        boundary = a_boundaries[x]*1000
        if(height <= boundary):
            layer = x
            break

    return layer


if __name__ == "__main__":
    startUp = True

    height = menuInput(startUp = True)

    layer = layerDetermination(height)

    print("The height selected was", str(height), "meters. This height indicates the", a_boundaries_names[layer], ".")


