allUnfiltered  =[76082, 97286, 85775]
allPrecip =  [1274, 1688, 1316]

allUnfiltered  =[68747, 70164, 69804]
allPrecip =  [475, 133, 157]


def percentIncorporation(allPrecip, allUnfiltered):	
	precipMean = sum(allPrecip)/len(allPrecip)
	unfilteredMean = sum(allUnfiltered)/ len(allUnfiltered) 
	num = precipMean * 65/45
	denom = unfilteredMean * 32.5
	answer = (num / denom) * 100
	return answer
percentIncorporation = percentIncorporation(allPrecip, allUnfiltered)


	
allUnfilteredRxn1  =[68747, 70164, 69804]
allPrecipRxn1 =  [475, 133, 157]

allErrorUnfilteredRxn1 = [.76, .76, .76]
allErrorPrecipRxn1 = [9.18, 17.34, 15.96]


	
allUnfiltered  =[76082, 97286, 85775]
allPrecip =  [1274, 1688, 1316]

allErrorUnfiltered = [0.73, .64, .68]
allErrorPrecip = [5.8, 4.87, 5.51]

rxn1stdevP = 190.9031168
rxn1stdevU = 736.5163497

rxn2stdevP = 227.8683831
rxn2stdevU = 10614.98144




import math
def errorIncorp(errorPrecip, errorUnfiltered, cpmPrecip, cpmUnfiltered):
	precipMean,unfilteredMean  = sum(allPrecip)/len(allPrecip), sum(allUnfiltered)/ len(allUnfiltered)
	prelim = float(precipMean) / float(unfilteredMean)
	answer = prelim * math.sqrt(( errorPrecip/ precipMean)**2 + ( errorUnfiltered/ unfilteredMean)**2)
	return answer	
errorIncorp = errorIncorp(rxn1stdevP, rxn1stdevU, allPrecip, allUnfiltered)


def errorIncorp(errorPrecip, errorUnfiltered, cpmPrecip, cpmUnfiltered):
    #errorPrecip, errorUnfiltered = sum(allErrorPrecip) / len(allErrorPrecip), sum(allErrorUnfiltered) / len(allErrorUnfiltered)
    prelim = errorPrecip / errorUnfiltered
    precipMean,unfilteredMean  = sum(allPrecip)/len(allPrecip), sum(allUnfiltered)/ len(allUnfiltered)
    answer = prelim * math.sqrt(( errorPrecip/ precipMean)**2 + ( errorUnfiltered/ unfilteredMean)**2)
    return answer	
errorIncorp = errorIncorp(rxn2stdevP, rxn2stdevU, allPrecip, allUnfiltered)