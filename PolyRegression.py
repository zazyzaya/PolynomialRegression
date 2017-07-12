# # # # # # # # # # # # # #
#  Polynomial Regression  #
#      GA Practice        #
#      Isaiah King        #
# # # # # # # # # # # # # #

import numpy
import pandas
import matplotlib.pyplot as plt
#import matplotlib.polynomial.polynomial as polynom
from random import *
import math

POP_SIZE = 1000
NUM_POINTS = 500
NUMCHILDREN = 4
MUTANTCOEF = 3
POLY_RANK = math.floor(random() * 10)
GENERATIONS = 10
RANDOM_ERROR = .3
INCRIMENTS = 100

def initializePoints (numPoints, polyRank):
    global NUM_POINTS
    
    pointArray = [[],[]]
    coefArray = []

    for i in range(polyRank + 1):
        c = 1
        if i == 0:
            c = 0
        randNum = random() + c*1
        if (random() > .5):
            randNum *= -1

        coefArray.append(randNum)

    for point in range(NUM_POINTS):
        x = random()
        pointArray[0].append(x)
        pointArray[1].append(plotOrigPoints(coefArray, x))

    return pointArray
    
def plotOrigPoints (coefArray, x):

    y = 0
    for power in range(len(coefArray)):
        y += ((x * coefArray[power]) ** power)

    jitter = RANDOM_ERROR * random()
    if random() > .5:
        jitter *= -1

    y += jitter
    return y

def plotPoints (coefArray, x):

    x = float(x)
    y = 0
    for power in range(len(coefArray)):
        try:
            y += ((x * coefArray[power]) ** power) 
        except (OverflowError):
            y += 99
    return y

def generatePopulation():
    popArray = []
    for i in range(POP_SIZE):
        numCoefs = math.floor((random() * 24)) + 2   # maximum of 25 degree polynomial
        gene = []
        for coef in range(numCoefs):
            notZero = 1
            if coef == 0:
                notZero = 0
            nucleotide = random() + 1*notZero
            if random() > .5:
                nucleotide *= -1
            gene.append(nucleotide)
        
        popArray.append([gene, 999])
    
    return popArray

def testFitness (gene, pointArray):
    
    fitness = 0
    for i in range(INCRIMENTS):
        x = i/INCRIMENTS
        y = plotPoints(gene[0], x)
        pArrayx = min(pointArray[0], key=lambda lam:abs(lam-x))     # finds closest value to x
        pArrayy = pointArray[1][pointArray[0].index(pArrayx)]          # Finds matching y value
        # d = distance(x, y, pArrayx, pArrayy)
        
        if (pArrayy - y) < 0:
            try:
                fitness += 99
            except (OverflowError):
                fitness = 99999
        else:
            fitness += abs(pArrayy - y)
    
    #fitness += abs(polyIntegral(gene[0]))/4
    return fitness

# Finds area of shape made by polynomial
def polyIntegral(coefs, lowBound=0, upBound=1):
    integralCoefs = coefs

    y1 = plotPoints(coefs, lowBound)
    y2 = plotPoints(coefs, upBound)

    a = distance(lowBound, y1, upBound, y2)
    b = upBound - lowBound

    # Find area benieth curve using trapazoid method
    topIntegral = ((a + b)/2) * max(y1, y2)

    for power in range(len(coefs)):
        integralCoefs.append(coefs[power]/(power+1))

    integralCoefs.insert(0,0)
    bottomIntegral = plotPoints(integralCoefs, upBound) - plotPoints(integralCoefs, lowBound)

    return topIntegral - bottomIntegral

def killWeaklings(pop):
    #print(len(pop)/NUMCHILDREN)
    pop = sorted(pop,key=lambda l:l[1])   # Sorts by lowest to highest rank
    pop = pop[:int(len(pop)/NUMCHILDREN)]                   # Pulls out strongest 4th
    print(pop[0][1])
    return pop

def distance(x1, y1, x2, y2):
    try:
        return ((x1 - x2) ** 2 + (y1 - y2)** 2) ** (1/2)
    except (OverflowError):
        return (999)

def giveBirth(population):
    print(len(population))
    shuffle(population)
    i = 0
    retPopulation = []

    while(i < len(population)):
        male = population[i][0]
        maleLen = len(male)
        female = population[i+1][0]
        femaleLen = len(female)

        if (maleLen > femaleLen):
            placeholder = male
            male = female
            female = placeholder
            maleLen = len(male)
            femaleLen = len(female)

        for j in range(NUMCHILDREN * 2):
            child = []
            for nucleotide in range(randint(maleLen,femaleLen)):   # Can swap up to 10 times
                takeFromMale = (1 == randint(0,1))

                if takeFromMale == True:
                    try:
                        child.append(male[nucleotide])
                    except IndexError:
                        break
                else:
                    try:
                        child.append(female[nucleotide])
                    except IndexError:
                        break
            
            if randint(0, MUTANTCOEF) == MUTANTCOEF:
                child = mutate(child)

            retPopulation.append([child, 999])
        i += 2
    
    print(len(retPopulation))
    return retPopulation

def mutate(dna):
    
    howManyMutations = randint(1, 2)
    for i in range(howManyMutations):
        whichMutation = randint(0, 2)
        
        # Replacement
        if whichMutation == 0:
            nucleotide = int(randint(0, len(dna)-1))
            dna[nucleotide] += (random() - .5) * 4

        # Addition
        elif whichMutation == 1:
            dna.append((random() - .5) * 2)

        # Subtraction
        elif whichMutation == 2:
            dna[int(randint(0, len(dna)-1))] = 0

    return dna

def main():
    pointArray = initializePoints(NUM_POINTS, POLY_RANK)
    population = generatePopulation()
    #topDogs = []
    best = [[1], 99999]

    for i in range(GENERATIONS):
        for gene in population:
            gene[1] = testFitness(gene, pointArray)
        
        currentBest = sorted(population,key=lambda l:l[1])[0]
        
        if currentBest[1] < best[1]:
            best = currentBest
        else:
            population.pop(POP_SIZE - 1)
            population.append(best)

        if best[1] <= 7:
            break
        population = killWeaklings(population)
        #topDogs.append(best)
        population = giveBirth(population)
        print(i)

    xArray = numpy.arange(0, 1, .01)
    yArray = []
    #for topDog in topDogs:
    for xVal in xArray:
        yArray.append(plotPoints(best[0], xVal))
    plt.plot(xArray, yArray, '-')
    #    yArray = []
    
    plt.plot(pointArray[0], pointArray[1], '.')
    plt.ylim(ymin = -5)
    plt.ylim(ymax = 5)
    plt.show()
       
main()
# while (True):
#    points = initializePoints(NUM_POINTS, POLY_RANK)
#    plt.plot(points[0], points[1], '.')
#    plt.ylim(ymin=-1)
#    plt.ylim(ymax=1)
#    plt.show()