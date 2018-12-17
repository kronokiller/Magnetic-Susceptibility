from statistics import stdev
from numpy import linspace
from math import log10, floor

# Measurement uncertainties
dmassMeasurements = .002
dheightMeasurements = 1
dcurrentMeasurements = .02

# Unmeasured values
cuvetteArea = 100
stripLength = 10
g = 9810
permeabilityOfFreeSpace = 1.25663706

#data collected
sampleMasses = [52.76, 35.23, 35.23, 41.164, 44.519, 27.813, 28.904, 32.149, 33.788, 36.472, 31.344, 30.504, 32.475, 30.831, 33.097, 31.743, 32.415]
sampleHeights = [32, 31, 31, 31, 32, 31, 29, 34, 32, 37, 32, 37, 32, 31, 33, 31]
cuvetteMasses = [27.553, 27.549, 27.556, 27.550]
calibrationMasses = [-.865, -.823, -.777, -.731, -.688, -.644, -.602, -.560, -.514, -.468, -.428, -.379, -.337, -.295, -.252, -.205, -.160, -.118, -.078, -.031, .01, .056, .099, .139, .185, .231, .274, .317, .360, .402, .441, .490, .534, .580, .621, .663, .707, .751, .794, .840, .883, .875, .830, .786, .742, .698, .654, .610, .568, .525, .480, .437, .395, .352, .310, .267, .221, .178, .137, .093, .050, .007, -.08, -.079, -.119, -.164, -.209, -.231, -.295, -.338, -.383, -.425, -.469, -.510, -.555, -.595, -.640, -.684, -.729, -.772, -.817, -.856]
# Second set of currents reversed to account for rotating the probe 180 degrees
calibrationCurrents = list(linspace(-2, 2, num = 41)) + list(linspace(2, -2, num = 41))
massChangesSampleTop = [.007, -.011, -.136, .071, -.203, .088, -.200, -2.499, -5.240, -.348, -.306, -.078, -.037, -1.204, -.860, -.891]
massChangesSampleBottom = [-.008, .012, .141, -.066, .175, -.087, .2, 2.569, 5.280, .313, .277, .078, .040, 1.211, .803, .767]

# Calculation of mass of empty cuvettes by taking mean and standard deviation
cuvetteMass = sum(cuvetteMasses) / len(cuvetteMasses)
dcuvetteMass = stdev(cuvetteMasses)

# Calculation of tared sample masses with uncertainties added in quadrature
dtaredMasses = ((dmassMeasurements) ** (2) + (dcuvetteMass) ** (2)) ** (1 / 2)
taredMasses = []
for i in range(16):
    taredMasses.append(sampleMasses[i] - cuvetteMass)

# Remove currents less than 0.5 A to prevent precision loss
popNum = 0
for i in range(len(calibrationCurrents)):
    if calibrationCurrents[i - popNum] < 0.49 and calibrationCurrents[i - popNum] > -0.49:
        calibrationCurrents.pop(i - popNum)
        calibrationMasses.pop(i - popNum)
        popNum += 1

# Calculation of each value of H
Hvalues = []
for i in range(len(calibrationMasses)):
    Hvalues.append(calibrationMasses[i] * g / permeabilityOfFreeSpace / stripLength / calibrationCurrents[i])

# Calculation of H where H is the mean dH is the standard deviation
H = sum(Hvalues) / len(Hvalues)
dH = stdev(Hvalues)

# Calculation of average mass changes the dmassChanges is the uncertainty in mass measurements for 2 measurements added in quadrature
massChanges = []
dmassChanges = 0.003
for i in range(len(massChangesSampleTop)):
    massChanges.append((-massChangesSampleTop[i] + massChangesSampleBottom[i]) / 2)

# Calculation of susceptibilities with the uncertaintiees in mass and H added in quadrature for each value of the susceptibility
susceptibilities = []
dsusceptibilities = []
for i in range(16):
    susceptibilities.append(2 * permeabilityOfFreeSpace * massChanges[i] * g / cuvetteArea / H ** 2)
    dsusceptibilities.append(((dmassChanges * susceptibilities[i] / massChanges[i]) ** 2 + (dH * 2  * susceptibilities[i] / H) ** 2) ** (1 / 2))
    print(dH/H, dmassChanges/massChanges[i])
    precision = floor(log10(dsusceptibilities[i]))
    print('(' + str(int(round(susceptibilities[i] / 10 ** precision))), '+/-', str(int(round(dsusceptibilities[i] / 10 ** precision))) + ')', 'E' + str(precision))