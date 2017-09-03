#!/usr/bin/env python 
import sys
import math
from Em import estMax

## Task 1

def parseDataFile(dataFile):
    runCount = 1
    data = []
    dataLine = dataFile.readline()
    while dataLine:
        if dataLine == "\n":  # Skip blank lines and assume new run
            dataLine = dataFile.readline()
            runCount += 1
            continue

        x = int(dataLine[1])
        y = int(dataLine[3])
        emissionIndex = int(dataLine[5:]) + 1  # Add 1 to get 0 -> 2 index

        data.append([runCount, [x, y], emissionIndex])

        dataLine = dataFile.readline()

    return data


def createGrid(value = 0):
    grid = []
    for i in range(0, 4):
        grid.append([value, value, value, value])

    return grid


def calculateInitialMles(data):
    runCount = 0
    grid = createGrid()

    for state in data:
        if state[0] == runCount:
            continue

        grid[state[1][0]][state[1][1]] += 1 # Record initial grid position
        runCount += 1

    mles = [[(count/runCount) for count in y] for y in grid]

    # format and print results
    print('\nInitials:')
    for x in range(4):
        for y in range(4):
            print('p(h1 = (', x, ',', y, '))  =  ', round(mles[x][y], 6))


def calculateEmissionMles(data):

    grid = createGrid()

    for state in data:
        x = state[1][0]
        y = state[1][1]
        emissionIndex = state[2]
        if grid[x][y] == 0:
            grid[x][y] = [0,0,0]

        grid[x][y][emissionIndex] += 1

    mles = grid[:]

    for xIndex, x in enumerate(mles):
        for yIndex, y in enumerate(x):
            if y == 0:
                mles[xIndex][yIndex] = '-'
                continue

            yTotalEmissions = sum(y)
            for emissionIndex, emissionCount in enumerate(y):
                mles[xIndex][yIndex][emissionIndex] = emissionCount/yTotalEmissions


    # format and print results
    print('\nEmissions:')
    for x in range(4):
        for y in range(4):
            for v in range(3):
                if v == 0: # Just for alignment of output
                    space = ''
                else:
                    space = ' '

                if mles[x][y] == '-':
                    print('p(vt = ', space, v-1,' | ht = (', x, ',', y, '))  =  ',  'Undefined')
                else:
                    print('p(vt = ', space, v-1,' | ht = (', x, ',', y, '))  =  ', round(mles[x][y][v], 6))


def calculateTransisionMles(data):
    grid = createGrid()
    runCount = 0

    for state in data:
        if state[0] != runCount:
            runCount += 1
            previousState = state
            continue

        x = state[1][0]
        y = state[1][1]
        previousX = previousState[1][0]
        previousY = previousState[1][1]
        if grid[previousX][previousY] == 0:
            grid[previousX][previousY] = [0,0,0,0,0] # Transition counts [LEFT, RIGHT, DOWN, UP, LOOP]

        if x == previousX - 1:
            grid[previousX][previousY][0] += 1
        elif x == previousX + 1:
            grid[previousX][previousY][1] += 1
        elif y == previousY - 1:
            grid[previousX][previousY][2] += 1
        elif y == previousY + 1:
            grid[previousX][previousY][3] += 1
        elif x == previousX and y == previousY:
            grid[previousX][previousY][4] += 1
        else:
            print('ERROR DATA WRONG')
        previousState = state

    mles = grid[:]


    # normalise counts
    for xIndex, x in enumerate(mles):
        for yIndex, y in enumerate(x):
            if y == 0:
                mles[xIndex][yIndex] = ['Undefined','Undefined','Undefined','Undefined','Undefined']
                continue

            yTotalTransitions = sum(y)
            for transisionIndex, transisionCount in enumerate(y):
                mles[xIndex][yIndex][transisionIndex] = round(transisionCount/yTotalTransitions, 6)

    # format and print results
    print('\nTransitions: ')
    for x in range(4):
        for y in range(4):

            #LEFT
            nsX = x-1
            nsY = y
            if nsX >= 0:
                print('p(ht+1 = (', nsX, ',', nsY, ') | ht = (', x, ',', y, '))  =  ', mles[x][y][0])

            #RIGHT
            nsX = x+1
            nsY = y
            if nsX < 4:
                print('p(ht+1 = (', nsX, ',', nsY, ') | ht = (', x, ',', y, '))  =  ', mles[x][y][1])

            #DOWN
            nsX = x
            nsY = y-1
            if nsY >= 0:
                print('p(ht+1 = (', nsX, ',', nsY, ') | ht = (', x, ',', y, '))  =  ', mles[x][y][2])

            #UP
            nsX = x
            nsY = y+1
            if nsY < 4:
                print('p(ht+1 = (', nsX, ',', nsY, ') | ht = (', x, ',', y, '))  =  ', mles[x][y][3])

            #LOOP
            print('p(ht+1 = (', x, ',', y, ') | ht = (', x, ',', y, '))  =  ', mles[x][y][4])

def task1(dataFilePath):
    dataFile = open(dataFilePath, 'r')
    data = parseDataFile(dataFile)

    calculateInitialMles(data)
    calculateEmissionMles(data)
    calculateTransisionMles(data)

## Task2-4

def parseEmissionData(dataFile):
    index = 0
    data = [[]]
    dataLine = dataFile.readline()
    while dataLine:
        if dataLine == "\n":  # Skip blank lines and assume new run
            dataLine = dataFile.readline()

            if dataLine:
                data.append([])
            index += 1
            continue

        emissionIndex = int(dataLine) + 1  # Add 1 to get 0 -> 2 index
        data[index].append(emissionIndex)

        dataLine = dataFile.readline()

    return data

def printModelParameters(initials, emissions, transitions, allowedTransisions = False):
    print('\nInitials:')
    cell = 0
    for y in range(4):
        for x in range(4):
            print('p(h1 = (', x, ',', y, '))  =  ', round(initials[cell], 6))
            cell += 1

    print('\nEmissions:')
    cell = 0
    for y in range(4):
        for x in range(4):
            for v in range(3):
                if v == 0: # Just for alignment of output
                    space = ''
                else:
                    space = ' '
                    
                print('p(vt = ', space, v-1,' | ht = (', x, ',', y, '))  =  ', round(emissions[cell][v], 6))
            cell += 1

    print('\nTransitions:')
    cell = 0
    for y in range(4):
        for x in range(4):
            nextState = 0
            for nsY in range(4):
                for nsX in range(4):
                    if allowedTransisions:
                        if allowedTransisions[cell][nextState]:
                            print('p(ht+1 = (', nsX, ',', nsY, ') | ht = (', x, ',', y, '))  =  ', round(transitions[cell][nextState], 6))
                    else:
                        print('p(ht+1 = (', nsX, ',', nsY, ') | ht = (', x, ',', y, '))  =  ', round(transitions[cell][nextState], 6))
                    nextState +=1
            cell += 1


def task2(dataFilePath):
    dataFile     = open(dataFilePath, 'r')
    observations = parseEmissionData(dataFile)

    Em = estMax(observations)
    Em.baumWelch()

    printModelParameters(Em.initials, Em.emissions, Em.transitions)


def task3(dataFilePath):
    dataFile     = open(dataFilePath, 'r')
    observations = parseEmissionData(dataFile)

    Em = estMax(observations)
    for i in range(10):
        print('\n\n\nRUN ', i+1)

        Em.randomise()
        Em.baumWelch()

        printModelParameters(Em.initials, Em.emissions, Em.transitions)


def task4(dataFilePath):
    dataFile     = open(dataFilePath, 'r')
    observations = parseEmissionData(dataFile)

    Em = estMax(observations)
    for i in range(10):
        print('\n\n\nRUN ', i + 1)

        Em.setConstraints()
        Em.randomise()
        Em.baumWelch()

        printModelParameters(Em.initials, Em.emissions, Em.transitions, Em.permittedTransisions)


if sys.argv[1] == '1':
    task1(sys.argv[2])

if sys.argv[1] == '2':
    task2(sys.argv[2])

if sys.argv[1] == '3':
    task3(sys.argv[2])

if sys.argv[1] == '4':
    task4(sys.argv[2])
