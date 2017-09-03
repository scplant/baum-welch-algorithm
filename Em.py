import random
import math

class estMax(object):

    def __init__(self, observations):
        self.numStates  = 16
        self.initials     = [1/self.numStates for s in range(0, self.numStates)]
        self.transitions  = [[1/self.numStates for t in range(0, self.numStates)] for s in range(0, self.numStates)] # all transition probabilities for all states
        self.emissions    = [[1/3 for e in range(0,3)] for s in range(0, self.numStates)]
        self.observations = observations

        self.permittedTransisions = [[True for t in range(0, self.numStates)] for s in range(0, self.numStates)]

    ##Helper functions


    # set the map of allowed transitions as specified for task 4
    def setConstraints(self):
        #     1      2    3      4         5      6      7     8        9     10      11    12        13     14     15     16
        self.permittedTransisions = [
            [True, True, False, False,   True, False, False, False,   False, False, False, False,   False, False, False, False],
            [True, True, False, False,   False, True, False, False,   False, False, False, False,   False, False, False, False],
            [False, False, True, True,   False, False, True, False,   False, False, False, False,   False, False, False, False],
            [False, False, True, True,   False, False, False, True,   False, False, False, False,   False, False, False, False],

            [True, False, False, False,   True, False, False, False,   True, False, False, False,    False, False, False, False],
            [False, True, False, False,   False, True, True, False,    False, True, False, False,    False, False, False, False],
            [False, False, True, False,   False, True, True, True,     False, False, True, False,    False, False, False, False],
            [False, False, False, True,   False, False, True, True,    False, False, False, False,   False, False, False, False],

            [False, False, False, False,  True, False, False, False,   True, False, False, False,    True, False, False, False],
            [False, False, False, False,  False, True, False, False,   False, True, False, False,    False, True, False, False],
            [False, False, False, False,  False, False, True, False,   False, False, True, True,     False, False, True, False],
            [False, False, False, False,  False, False, False, False,  False, False, True, True,     False, False, False, True],

            [False, False, False, False,  False, False, False, False,  True, False, False, False,    True, True, False, False],
            [False, False, False, False,  False, False, False, False,  False, True, False, False,    True, True, False, False],
            [False, False, False, False,  False, False, False, False,  False, False, True, False,    False, False, True, True],
            [False, False, False, False,  False, False, False, False,  False, False, False, True,    False, False, True, True]
            ]


    def randomise(self):
        self.initials     = [random.randint(0,1000) for s in range(0, self.numStates)]
        self.transitions  = [[random.randint(0,1000) for t in range(0, self.numStates)] for s in range(0, self.numStates)]
        self.emissions    = [[random.randint(0,1000) for e in range(0,3)] for s in range(0, self.numStates)]

        for cell in range(self.numStates): # limit to permitted transitions (All for task 3)
            for transition in range(self.numStates):
                if self.permittedTransisions[cell][transition] is False:
                    self.transitions[cell][transition] = 0

        self.normalise(self.initials)
        self.normalise(self.transitions)
        self.normalise(self.emissions)


    def createArray(self, states, columns):
        if states == 1:
            return [0 for c in range(0, columns)]

        return [[0 for c in range(0, columns)] for s in range(0, states)]


    # deep copy any list recursively
    def deepCopy(self, A):
        rt = []
        for elem in A:
            if isinstance(elem,list):
                rt.append(self.deepCopy(elem))
            else:
                rt.append(elem)
        return rt


    # normalise all values in lowest list nodes recursively
    def normalise(self, array):
        if isinstance(array[0], list):
            for elem in array:
                self.normalise(elem)
        else:
            total = 0.0
            for i in range(len(array)):
                total += array[i]
            for i in range(len(array)):
                array[i] = array[i]/total

    ## Algorithm
    def baumWelch(self):
        deltaLogLikelihood = 1
        newLogLikelihood   = 0

        while deltaLogLikelihood >= 0.01:
            oldInitials    = self.deepCopy(self.initials)
            oldTransitions = self.deepCopy(self.transitions)
            oldEmissions   = self.deepCopy(self.emissions)

            self.initials    = [0 for s in range(0, self.numStates)]
            self.transitions = [[0 for t in range(0, self.numStates)] for s in range(0, self.numStates)]
            self.emissions   = [[0 for e in range(0,3)] for s in range(0, self.numStates)]
            probOfObs = 1

            for observation in self.observations: # sum all observation sequences
                alpha = self.forward(observation, oldInitials, oldTransitions, oldEmissions)
                beta  = self.backward(observation, oldTransitions, oldEmissions)
                gamma, Xij = self.expectation(observation, oldTransitions, oldEmissions, alpha, beta)
		
                self.maximisation(observation, gamma, Xij)
                probOfObs *= sum([alpha[s][len(observation) - 1] for s in range(self.numStates)])

            self.normalise(self.initials)
            self.normalise(self.transitions)
            self.normalise(self.emissions)

            oldLogLikelihood = newLogLikelihood
            newLogLikelihood = math.log(probOfObs)
            deltaLogLikelihood = abs(newLogLikelihood - oldLogLikelihood)
            print('Log Likelihood: ', newLogLikelihood)


    # alpha probabilities P(o1, ..., ot| Theta)
    def forward(self, observations, initials, transisions, emissions):
        forwardProbs = self.createArray(self.numStates, len(observations))

        for s in range(self.numStates): # for t=0
            forwardProbs[s][0] = initials[s] * emissions[s][observations[0]]

        for t in range(1, len(observations)): # loop for all ts (observations)
            observation = observations[t]

            for s in range(self.numStates): # loop for all states
                sumPrevious = sum([forwardProbs[ps][t - 1] * transisions[ps][s] for ps in range(0, self.numStates)]) # sum total probability of reaching this state
                forwardProbs[s][t] = emissions[s][observation] * sumPrevious
                
        return forwardProbs


    # beta probabilities P(o(t+1), ..., oT | Theta)
    def backward(self, observations, transisions, emissions):
        backwardProbs = self.createArray(self.numStates, len(observations))

        for s in range(self.numStates): # for t = T
            backwardProbs[s][len(observations) - 1] = 1.0 # set final state to 1 (assume this state is reached)

        for t in range(len(observations) - 1, 0, -1): # loop from end to start of sequence
            observation = observations[t]

            for s in range(self.numStates): # loop for all states
                backwardProbs[s][t-1] = sum([transisions[s][ns] * emissions[ns][observation] * backwardProbs[ns][t] for ns in range(self.numStates)]) #ns is next state

        return backwardProbs


    def expectation(self, observations, transisions, emissions, alpha, beta):

        totalProbabilityObservation = sum([alpha[s][len(observations) - 1] for s in range(self.numStates)]) # P(X| Theta)

        gamma = self.createArray(self.numStates, len(observations))
        Xij   = self.createArray(self.numStates, self.numStates)
        for s in range(self.numStates):
            for ns in range(self.numStates):
                Xij[s][ns] = self.createArray(1, len(observations) - 1)

        for t in range(len(observations)):
            for s in range(self.numStates):
                gamma[s][t] = (alpha[s][t] * beta[s][t]) / totalProbabilityObservation

                if t == (len(observations) - 1): # compute  Xij up until T-1
                    continue

                for ns in range(self.numStates):
                    Xij[s][ns][t] = (alpha[s][t] * transisions[s][ns] * emissions[ns][observations[t+1]] * beta[ns][t+1])/ totalProbabilityObservation

        return (gamma, Xij)

    def maximisation(self, observations, gamma, Xij):

        for s in range(self.numStates):
            self.initials[s] += gamma[s][0] #calculate new initials

        for s in range(self.numStates): #calculate new transisions
            for ns in range(self.numStates):
                probT = 0.0
                for t in range(len(observations) - 1):
                    probT += Xij[s][ns][t]
                self.transitions[s][ns] += probT

        for s in range(self.numStates): #calculate new emissions
            for e in range(3): #loop over emissions
                probE = 0.0
                for t in range(len(observations)):
                    if observations[t] == e:
                        probE += gamma[s][t]
                self.emissions[s][e] += probE
