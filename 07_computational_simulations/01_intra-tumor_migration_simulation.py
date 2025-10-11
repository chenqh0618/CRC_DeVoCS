import sys
import os
os.chdir("/public/data/cqh_project/crc/simulation/inter_tumor_migration/")


#import packages
import numpy as np, math, time, random, multiprocessing,sys
from collections import Counter
import dill, copy, logging
import matplotlib ; matplotlib.use('Agg')
import seaborn as sns; sns.set_style('ticks')
import matplotlib.pyplot as plt
logging.basicConfig(format='%(levelname)s: %(processName)s, %(message)s', level=logging.DEBUG)

#timer
try:  # python 3.3 or later
    timer = time.perf_counter #(includes time during sleep)
    # timer = time.process_time #(not includes time during sleep)   
except AttributeError: # before python 3.3
    timer = time.clock if sys.platform[:3] == 'win' else time.time
    
def totalTime(reps, func, *pargs, **kargs):
    '''Total time to run func() reps times. Return (total time, last result)'''
    repslist = list(range(reps))
    start = timer()
    for i in repslist:
        ret = func(*pargs, **kargs)
    elapsed = timer() - start
    return (elapsed, ret)
    
def vonNeumannNeighbor(coord):
    """
    4 adjacent sites in 2D; 6 in 3D
    :: coord ::   tuple -- coordinates, like (x,y) or (x,y,z)
    :: return ::  set   -- coordinates of neiboring positions
    """
    if len(coord) == 2: # 2D
        x, y = coord
        neighbor = [(x-1, y), (x, y-1), (x, y+1), (x+1, y)]
    else:               # 3D
        x, y, z = coord
        neighbor = [(x-1, y, z), (x, y-1, z), (x, y+1, z), (x, y, z-1), (x, y, z+1), (x+1, y, z)]
    return neighbor
    
def MooreNeighbor(coord):
    """
    8 adjacent sites in 2D; 26 in 3D
    :: coord ::   tuple -- coordinates, like (x,y) or (x,y,z)
    :: return ::  set   -- coordinates of neiboring positions
    """
    if len(coord) == 2:  # 2D
        x, y = coord
        neighbor = {(x-1, y-1), (x-1, y), (x-1, y+1), (x, y-1), (x, y+1), (x+1, y-1), (x+1, y), (x+1, y+1)}
    else:                # 3D
        x, y, z = coord
        neighbor = {(x-1, y-1, z-1), (x-1, y-1, z), (x-1, y-1, z+1),
                    (x-1, y, z-1), (x-1, y, z), (x-1, y, z+1), 
                    (x-1, y+1, z-1), (x-1, y+1, z), (x-1, y+1, z+1), 
                    (x, y-1, z-1), (x, y-1, z), (x, y-1, z+1), 
                    (x, y, z-1), (x, y, z+1), (x, y+1, z-1), 
                    (x, y+1, z), (x, y+1, z+1), (x+1, y-1, z-1), 
                    (x+1, y-1, z), (x+1, y-1, z+1), (x+1, y, z-1),
                    (x+1, y, z), (x+1, y, z+1), (x+1, y+1, z-1),
                    (x+1, y+1, z), (x+1, y+1, z+1) 
                    }
    return neighbor
     
def localNeighbor(coord, r):
    """
    search the local neighbour sites of coord within a spherical area of radius r (including r) in the 2D or 3D cubic lattice.
    """
    iterator = range(-r, r+1)
    if len(coord) == 2: # 2D
        x, y = coord
        neighbor = {(x+i, y+j) for i in iterator for j in iterator if i**2 + j**2 <= r**2}
    else:               # 3D
        x, y, z = coord
        neighbor = {(x+i, y+j, z+k) for i in iterator for j in iterator for k in iterator if i**2 + j**2 + k**2 <= r**2}
    return neighbor      

def centralPlane(center, coordList3D):
    '''
    Get the coordinates in the cross section (through the center with z=z0) according to a list of 3D coordinates (coordList3D)
    '''
    z0 = center[2]
    crossSection = {coord for coord in coordList3D if coord[2] == z0}
    return crossSection
    
def traceLineage(mlineage, recentMut): ## recentMut is a string-format integer, e.g. '11', '11,12'
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell. 
    We use the tree data structure to record the mutations. For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is [100, 56, 12, 1], which is the mutation lineage of the cell
    ::mlineage::  list   -- a tree data structure that could be used to recover the mutational lineage given the most recent mutation id of a lineage, the element is a string-format integer, e.g. '1' or '2,3'
    ::recentMut:: string -- the mutation ID of the most recently occurred mutation in the cell, BCF: if the recentMut is '0', the function will return a null list []
    ::return::    list   -- list of integer like [100, 56, 12, 1], BCF: the element is integer type rather than string type, not containing the initial mutation id (0)
    """
    recent_muts = recentMut.split(',')  # it is possible that multiple mutations occur during in a cell division, e.g. "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0]          # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace += recent_muts
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts] 
        first_mut = recent_muts[0]
    return trace  
    
def multinomial(bProb, dProb, mProb):
    '''
    Revise multinomial verdict to the np.random.multinomial error pval[:-1] > 1
    '''
    maxProb = max(bProb, dProb, mProb)
    if maxProb < 1:
        isB = np.random.binomial(1, bProb)
        if isB: isD, isM, isPass = 0, 0, 0
        else:
            isD = np.random.binomial(1, dProb)
            if isD: isM, isPass = 0, 0
            else:
                isM =  np.random.binomial(1, mProb)
                isPass = 0 if isM else 1
    elif maxProb == 1: 
        if bProb == 1: isB, isD, isM, isPass = 1, 0, 0, 0
        elif dProb == 1: isB, isD, isM, isPass = 0, 1, 0, 0
        else: isB, isD, isM, isPass = 0, 0, 1, 0
    else: raise ValueError('maximal probability is not supposed to greater than 1')    
    return isPass, isB, isD, isM
    
class Cell(object):
    def __init__(self, birthR, deathR, recentMut):
        '''
        :: birthR :: float -- birth rate of a cell for continuous model.
        :: deathR :: float -- death rate of a cell for continuous model.
        :: recentMut :: string -- the mutation ID of the most recently occurred mutation in the cell. It can have multiple mutations for one divisiion, e.g. '2,3'
        '''
        self.b = birthR
        self.d = deathR
        self.recentMut = recentMut
        
class SpaceCell(Cell):
    def __init__(self, birthR, deathR, migrationR, recentMut, position):
        '''
        :: birthR :: float -- intrinsic birth rate of a cell which only determined by genetic genotype.
        :: deathR :: float -- intrinsic death rate of a cell which only determined by genetic genotype.
        :: recentMut :: string -- the mutation ID of the most recently occurred mutation in the cell. It can have multiple mutations for one divisiion, e.g. '2,3'
        :: envB, envD, envM :: float  -- the eventual birthR, deathR, migrationR affected by both genetic genotype and environmental factors. The default values are the same as intrinsic rates, i.e. not considering microenvironment
        '''
        Cell.__init__(self, birthR, deathR, recentMut)
        self.m = migrationR
        self.pos = position
        self.envB = birthR
        self.envD = deathR
        self.envM = migrationR
    def sumEnvR(self):
        return self.envB + self.envD + self.envM
    def sumR(self):
        return self.b + self.d + self.m
    def maxEnvR(self):
        return max(self.envB, self.envD, self.envM)
        
class Population(dict):
    def __init__(self, initCells, initT, u = None, 
                 posR = 0, posS = 0, negR = 0, negS = 0, 
                 mlineage = None, initMutID = 0, 
                 neighborFunc = MooreNeighbor, migrationPOS=localNeighbor):
        '''
        initial population has the following attributes:
        :: N     -- size of the population, i.e. the number of cells in the cell population
        :: countID  -- the total number of cells having occurred, include the dead cells. We mark each cell by a unique id (the keys, from 1 to countID).
        :: mutID    -- the total number of mutations having occurred, include the extinct mutations because of cell death. We mark each mutation by a unique id (from initMutID to mutID).
        :: curT     -- the time the poulation evolves
        :: u        -- the mutation rate per exome per cell division. u is considered to be constant unless we including the mutator mutation. if u is None (default) or 0, we don't consider mutation, just cell divisions
        :: posR, negR -- the probability that a CDS mutation is positive mutation, negative mutation respectively. if 0 (default), we don't consider it.
        :: posS, negS -- nonnegative number, selction coefficient for a positive, negative mutation.
        :: mlineage  -- a list that could be used to recover the mutational lineage given the most recent mutation id of a lineage, the element is a string-format integer, e.g. '1' or '2,3'
        :: neighborFunc -- function, determine the neighborhood of a position. The default is MooreNeighbor. An alternative option is vonNeumannNeighbor. Or we can define other neighborhood.
        BCF: we only use string-format integer to record mutations in mlineage. In posmutL, negmutL, or the return of traceLineage function, we all directly use interger to record mutations. Because the former will have multiple mutations, the latter won't.
        '''
        initSize = len(initCells)
        self.N = initSize
        self.countID = initSize
        self.mutID = initMutID  # the initMutID represents the number of mutations in the initial population, default 0.
        self.curT = initT
        self.u = u
        self.posR, self.negR = posR, negR
        self.posS, self.negS = posS, negS
        self.mlineage = ['0'] if mlineage is None else mlineage    # mlineage is mutable object (list), can't use self.mlineage = mlineage and self initial default mlineage = ['0']
        self.posmutL = [] 
        self.negmutL = [] 
        self.neighborFunc = neighborFunc
        self.migrationPOS=migrationPOS
        for cell in initCells:
            self[cell.pos] = cell
    def size(self):
        return len(self)
    def update(self):
        ''' update the attributes of all the cells in the cell population, return None '''
        for pos, cell in self.items():
            vacantPos, occupyPos = self.divideNeighborPos(pos)
            vacantNum = len(vacantPos)
            self.updateCellAttr(cell, vacantNum)
    def maxEnvR(self):
        ''' update envR of all the cells, and return the maximal envR after update '''
        maxER = 0
        for pos, cell in self.items():
            vacantPos, occupyPos = self.divideNeighborPos(pos)
            vacantNum = len(vacantPos)
            self.updateCellAttr(cell, vacantNum)
            maxER = max(maxER, cell.maxEnvR())
        return maxER
    def totEnvR(self):         
        ''' update envR of all the cells, and return the totEnvR after update '''
        totER = 0
        for pos, cell in self.items():
            vacantPos, occupyPos = self.divideNeighborPos(pos)
            vacantNum = len(vacantPos)
            self.updateCellAttr(cell, vacantNum)
            totER += (cell.envB + cell.envD + cell.envM)  
        return totER

    def posEffect(self, cell):
        '''
        change a cell's attributes according to positive mutation effect. 
        Here we assume postive mutation can only increase birth rate. More effect can be expanded.
        '''
        cell.b *= (1 + self.posS)
    def negEffect(self, cell):
        '''
        change a cell's attributes according to negative mutation effect. 
        Here we assume negative mutation can only increase death rate. More effect can be expanded.
        '''
        cell.d *= (1 + self.negS)
    def IAM(self): 
        '''
        Caculate the spectrum of allele according to the infinite allele model. Here a recentMut of a cell represents an allele.
        :: return :: dict   -- key is a unique allele represented by recentMut, value is the frequency of this allele. e.g. Counter({'1': 2, '2': 1, '3,4': 1}): there are 3 alleles with frequencies 2, 1, 1 respectively.
        '''
        allele_list = [cell.recentMut for cell in self.values()] # the recentMut of a cell reprsent an allele
        cellNum = len(allele_list)
        allele_count = Counter(allele_list)
        return allele_count
    def ISM(self):
        '''
        Caculate the spectrum of site according to the infinite site model. Here we get all the variant sites according to recentMut of each cell and then count the number of each variant site.
        :: return :: dict   -- key is a unique variatn site (str-format, like '1'), value is the frequency this variant site. e.g. Counter({'1': 2, '2': 1, '4': 1}): there are total 3 varian sites.
        '''
        allele_list = [cell.recentMut for cell in self.values()]
        sitePool = []
        for allele in allele_list:  # an allele contains at least one mutations, count the number of mutations
            lineage = traceLineage(self.mlineage, allele)
            sitePool.extend(lineage)
        site_count = Counter(sitePool)
        return site_count
    def mutCompound(self, cell):
        '''
        Obtain the number of positive, negative, neutral mutations respectively of a cell.
        :: return :: tuple -- the number of positive, negative, neutral mutations
        '''
        clineage = traceLineage(self.mlineage, cell.recentMut)
        totalMut = len(clineage)
        posL, negL = [], []
        for mut in clineage:
            # mut = str(mut)
            if mut in self.posmutL:
                posL.append(mut)
            elif mut in self.negmutL:
                negL.append(mut)
        posNum, negNum = len(posL), len(negL)
        return posNum, negNum, totalMut - posNum - negNum
    def divideNeighborPos(self, pos):  
        '''
        divide the neighborPos into vacant and occupied, then return.
        :: return :: tuple of two sets, the first is the vacant positions, the second is the occupied positions
        '''
        neighborPos = self.neighborFunc(pos)
        vacantPos = {pos for pos in neighborPos if pos not in self.keys()}
        occupyPos = neighborPos - vacantPos
        return vacantPos, occupyPos
    def updateCellAttr(self, cell, vacantNum): 
        '''
        update a cell's attributes according to its microenvironment. 
        # TODO: Here we only consider that the number of of vacant grids affect the birthR, more on death rate or migration rate need to be added
        '''
        if vacantNum: # have at least one vacant grids
            cell.envB = cell.b
        else:   # fullly surrounded by neiboring cells
            cell.envB = 0
    def updateLocalEnvR(self, posL): 
        '''
        update envR of cells in the local positions, and return some important rate.
        :: return :: tuple of rates -- totEnvRdelta, maxEnvB, maxEnvD, maxEnvM, maxEnvR of these chosen cells.
        '''
        totEnvRdelta, maxEnvB, maxEnvD, maxEnvM = 0, 0, 0, 0
        for pos in posL:
            if pos in self.keys():
                cell = self[pos]
                preEnvR = cell.sumEnvR()
                vacantPos, occupyPos = self.divideNeighborPos(pos)
                vacantNum = len(vacantPos)
                self.updateCellAttr(cell, vacantNum) # update attributes
                envB, envD, envM = cell.envB, cell. envD, cell.envM
                maxEnvB = envB if envB > maxEnvB else maxEnvB
                maxEnvD = envD if envD > maxEnvD else maxEnvD
                maxEnvM = envM if envM > maxEnvM else maxEnvM
                updateEnvR = envB + envD + envM
                totEnvRdelta += (updateEnvR - preEnvR)
        maxEnvR = max(maxEnvB, maxEnvD, maxEnvM)
        return totEnvRdelta, maxEnvB, maxEnvD, maxEnvM, maxEnvR
    def cellDivide(self, cell):
        '''
        A cell divides into two cells: one is the parental cell, the other is the offspring cell. The parental cell is just adjusted attributes from the input cell because of mutations. Thus, it's in the cell population after division. However, the second new cell is a totally new cell, it's not in the cell population. 
        Here, we only generate and return two cells assuming that the input cell can divide (In addition add the new mutations to the mlineage, posmutL, negmutL). The update of the cell population need to be carried out later (like adding the new cell to the population, update local rate)
        :: cell ::    instance of Cell class -- the input cell will divide into two cells
        :: return ::  tuple of two cells, the first is the parental cell, the other is the new cell
        '''
        # self.update()
        vacantPos, occupyPos = self.divideNeighborPos(cell.pos)
        print()
        targetPos = random.choice(list(vacantPos))
        b, d, m, recentMut = cell.b, cell.d, cell.m, cell.recentMut
        ## just let the attributes of newCell same as the parental cell first. We will change its attributes and the parent cell's attributes if there is new mutation to affect their attributes
        newCell = SpaceCell(b, d, m, recentMut, targetPos)
        if self.u:              # consider the mutation
            new_mut_num = np.random.poisson(lam = self.u * 2, size = None)  
            if new_mut_num == 0: 
            	pass  # the new mutation number is 0, the same as not considering mutations
            else:
                parent_mut_num, offspring_mut_num = np.random.multinomial(new_mut_num, [0.5, 0.5])
                ## process parent cell and update its attributes 
                parent_recentMut = ''
                for mut in range(parent_mut_num):  # parent cell has new mutations. if there is no new mutation (parent_mut_num == 0), just skip this for loop
                    self.mutID += 1
                    self.mlineage.append(recentMut)
                    parent_recentMut += ',%s' %self.mutID
                    posNum, negNum, neuNum = np.random.multinomial(1, [self.posR, self.negR, 1 - self.posR - self.negR])
                    if neuNum:   # neutral mutation
                        pass
                    elif posNum: # positive mutation
                        self.posEffect(cell)
                        self.posmutL.append(self.mutID)
                    else:        # negative mutation
                        self.negEffect(cell)
                        self.negmutL.append(self.mutID)
                if parent_recentMut: 
                	cell.recentMut = parent_recentMut[1:] ## if parent cell has new mutations, the first str of parent_recentMut is ','. remove it.
                ## process offspring cell and update its attributes
                offspring_recentMut = ''
                for mut in range(offspring_mut_num): # offspring cell has new mutations
                    self.mutID += 1
                    self.mlineage.append(recentMut)
                    offspring_recentMut += ',%s' %self.mutID
                    posNum, negNum, neuNum = np.random.multinomial(1, [self.posR, self.negR, 1 - self.posR - self.negR])
                    if neuNum:   # neutral mutation
                        pass
                    elif posNum: # positive mutation
                        self.posEffect(newCell)
                        self.posmutL.append(self.mutID)
                    else:        # negative mutation
                        self.negEffect(newCell)
                        self.negmutL.append(self.mutID)
                if offspring_recentMut: 
                	newCell.recentMut = offspring_recentMut[1:] ## if offspring cell has new mutations, the first str of parent_recentMut is ','. remove it.
        else:                   # don't consider mutation
            pass
        return cell, newCell
        
    def cellMigrate(self, cell):
        '''
        Carry out migration behavior of a input cell and return the positions in which a cell's attributes will need to be updated because of the migration.
        Here we assume the input cell will migrate (i.e. the migration verdict is not in the function). 
        :: return ::  set -- positions which need to be updated for the input cell's migration.
        '''
        pos = cell.pos
        neighborPos = self.neighborFunc(pos)
        migrationPOS=self.migrationPOS(pos,3)
        migPos = random.choice(list(migrationPOS))
        updatePos = neighborPos.union(self.neighborFunc(migPos).union({migPos,pos}))
        if migPos in self.keys(): # the migPos is occupied, exchange the two cell positions
            nativeCell = self[migPos]
            nativeCell.pos = pos
            self[pos] = nativeCell
            cell.pos = migPos
            self[migPos] = cell
        else:
            cell.pos = migPos
            migCell = self.pop(pos)
            self[migPos] = migCell
        return updatePos
        
    def Gillespie_stepForward(self, totER):
        '''
        Use Gillespie algorithm (more precisely, Gillespie Direct Method Algorithm) to simulate one step of the cell population's evolution.
        :: return :: float -- the total rate of all events.
        '''
        deltaT = np.random.exponential(1.0 / totER)
        self.curT += deltaT
        u2 = 1-np.random.random_sample() # random value from uniform distribution (0,1] 
        ru2, sumR = totER * u2, 0
        for pos in set(self.keys()): # Gillespie algorithm, each step, only one event happens (divide, die or migrate)
            cell = self[pos]
            b, d, m, recentMut, envB, envD, envM = cell.b, cell.d, cell.m, cell.recentMut, cell.envB, cell.envD, cell.envM
            if sumR < ru2 <= sumR + envB:  # this cell will divide
                cell, newCell = self.cellDivide(cell)
                self.countID += 1
                self.N += 1
                targetPos = newCell.pos
                self[targetPos] = newCell
                totER += newCell.sumEnvR()
                updatePos = self.divideNeighborPos(targetPos)[1].union({targetPos}) # after division, this positions need to be updated, which will include the position of parental cell
                totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
                totER += totEnvRdelta
                return totER  # break the for loop by return
            else: sumR += envB
            if sumR < ru2 <= sumR + envD:  # this cell will die
                self.N -= 1
                deadCell = self.pop(pos)
                totER -= deadCell.sumEnvR()
                updatePos = self.neighborFunc(pos)
                totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
                totER += totEnvRdelta
                return totER  # break the for loop by return
            else: sumR += envD
            if sumR < ru2 <= sumR + envM:  # this cell will migrate
                updatePos = self.cellMigrate(cell)  # return the positions which need to be updated for the cell's migration event
                totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
                totER += totEnvRdelta
                return totER # break the for loop by return
            else: sumR += envM 
    def rKMC_stepForward(self, maxEnvR):
        '''
        Use rejection kinetic Monte Carlo(rKMC) method to simulate one step of the cell population's evolution.
        :: return :: float -- the maximum rate of all events.
        '''   
        N = self.N  # N = self.size()
        deltaT = np.random.exponential(1.0 / N / maxEnvR)
        self.curT += deltaT
        pos = random.choice(list(self.keys()))  # randomly pick a cell from cell population
        cell = self[pos]
        b, d, m, recentMut, envB, envD, envM = cell.b, cell.d, cell.m, cell.recentMut, cell.envB, cell.envD, cell.envM
        bProb, dProb, mProb = envB / maxEnvR, envD / maxEnvR, envM / maxEnvR
        isPass, isB, isD, isM = multinomial(bProb, dProb, mProb)
        if isPass:   # No event is accepted, i.e. no birth or death
            return maxEnvR
        elif isB:    # the chosen cell will divide  
            cell, newCell = self.cellDivide(cell)
            self.countID += 1
            self.N += 1
            targetPos = newCell.pos
            self[targetPos] = newCell
            occupyPos = self.divideNeighborPos(targetPos)[1]  # the occupied pos of the neighborPos of the targetPos
            updatePos = occupyPos.union({targetPos,pos})
            #### use a suitable upper bound rate as maxEnvR. This rate can be not the exact maximal rate of all the cells, but it must be not less than the exact maximal rate. It's equivalent to the following method, but with lower accept probability.
            totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
            maxEnvR = max(maxEnvR, local_maxEnvR)
            #### use the exact maximal rate of all the cells as maxEnvR.
            # premaxEnvR = max(self[i].maxEnvR() for i in occupyPos) # caculate the maxEnvR before update to determine if the maxEnvR is from these cells in occupyPos.
            # totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
            # maxEnvR = max(local_maxEnvR, maxEnvR)
            # if maxEnvR > premaxEnvR:
                # return maxEnvR     
            # else:   # TODO: maybe many cells have the maximal rate at the same time.
                # return self.maxEnvR()  # updating the whole population takes lots of time
            return maxEnvR
        elif isD:    # the chosen cell will die
            self.N -= 1
            deadCell = self.pop(pos) 
            updatePos = self.neighborFunc(pos)
            totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
            maxEnvR = max(local_maxEnvR, maxEnvR)
            # maxEnvR = maxEnvR if maxEnvR > deadCell.maxEnvR() else self.maxEnvR() # Maybe the maxEnvR is in dead cell's rates, but maxEnvR can be greater than the exact maximum rate. Thus, this command can be deleted.
            return maxEnvR
        elif isM:    # the chosen cell will migrate
            updatePos = self.cellMigrate(cell)  # return the positions which need to be updated for the cell's migration event
            totEnvRdelta, local_maxEnvB, local_maxEnvD, local_maxEnvM, local_maxEnvR = self.updateLocalEnvR(updatePos)
            maxEnvR = max(local_maxEnvR, maxEnvR) # Maybe the maxEnvR extinct , but maxEnvR can be greater than the exact maximum rate.
            return maxEnvR
        else: raise ValueError('The list of multinomial probability is incorrect')
                
    def Nowak_stepForward(self):
        # TODO
        pass
        
def Gillespie(popul, maxSize):
    np.random.seed()
    random.seed()
    N = popul.size()
    assert maxSize > N, "The current populatin size should be less than maxSize."
    tL, nL = [popul.curT], [N]
    totEnvR = popul.totEnvR()
    while N <= maxSize:
        totEnvR = popul.Gillespie_stepForward(totEnvR)
        N = popul.N
        if N % 100 == 0:
            logging.debug('Gillespie size: %s' %N)
        assert N == popul.size()
        # assert abs(totEnvR - popul.totEnvR()) <= 0.01 # take too much time when size is big
        tL.append(popul.curT)
        nL.append(N)
        if N == 0: ## population extinction
            return None
    return popul, tL, nL

def rKMC(popul, maxSize):
    np.random.seed()
    random.seed()
    N = popul.size()
    assert maxSize > N, "The current populatin size should be less than maxSize."
    tL, nL = [popul.curT], [N]
    maxEnvR = popul.maxEnvR()
    while N <= maxSize:
        maxEnvR = popul.rKMC_stepForward(maxEnvR)
        N = popul.N
        if N % 100 == 0:
            logging.debug('rKMC size: %s' %N)
        assert N == popul.size()
        # tmaxEnvR = popul.maxEnvR()  # take too much time when size is big
        # assert maxEnvR == tmaxEnvR, '%s\t%s' %(maxEnvR, tmaxEnvR)
        tL.append(popul.curT)
        nL.append(N)
        if N == 0: ## population extinction
            return None
    return popul, tL, nL   
    
def Nowak(popul, maxSize):
    # TODO
    pass
 
def growthRate_compare(maxS, reps, proc = 5, outName = 'growthRate_compare.jpg'):    
    b, d, m = 0.693, 0, 0.1
    net = b - d
    initP = Population(initCells = [SpaceCell(b, d, m, '0', (0,0))], initT = 0, u = 0.1, posR = 0.01, negR = 0, posS = 0.1, negS = 0, mlineage = None, initMutID = 0)
    N0 = initP.size()
    argsL = [(initP, maxS)]*reps
    pools = multiprocessing.Pool(processes = proc)
    # nL = pools.starmap_async(Nowak, argsL)
    rL = pools.starmap_async(rKMC, argsL)
    gL = pools.starmap_async(Gillespie, argsL)
    pools.close(); pools.join()
    # gL, nL, rL = [res for res in gL.get() if res], [res for res in nL.get() if res],  [res for res in rL.get() if res] # delete the extinct populations
    gL, rL = [res for res in gL.get() if res], [res for res in rL.get() if res] # delete the extinct populations
    fig, ax = plt.subplots()
    for i, res in enumerate(gL):
        xs, ys = res[1:]
        ax.plot(xs, ys, color = 'red', label = 'Gillespie' if i == 0 else '')
    ax.legend()
    # for i, res in enumerate(nL):
        # xs, ys = res[1:]
        # ax.plot(xs, ys, color = 'green', label = 'Nowak' if i == 0 else '')
    ax.legend()
    for i, res in enumerate(rL):
        xs, ys = res[1:]
        ax.plot(xs, ys, color = 'blue', label = 'rKMC' if i == 0 else '')
    ax.legend()
    maxT = np.log(maxS/N0) / net
    xs = np.linspace(0, maxT, 101)
    ys = N0 * np.exp(net*xs)
    ax.plot(xs, ys, color = 'black', label = 'Exponential', lw = 3)
    ax.legend()
    plt.savefig(outName)
    plt.close('all')
 
def SFS(i, totMutNum, norm = False):
    """
    Equation to get the site frequency spectrum according neutral evolution
    """
    ic = 1.0 / i / (i+1)
    if norm:
        return ic
    else:
        return ic * totMutNum

def sfs_compare_freq(siteL, n, stepNum = 20, maxFreq = 1, DPI = 300, outName = None):
    '''plot site frequency spectrum difference between simulation and equation'''
    labelL = ['Equation SFS', 'Simulation SFS']
    mutNum, maxCount = len(siteL), int(round(n*maxFreq))
    xs = np.linspace(0, maxFreq, stepNum+1)
    xticks = (maxCount*xs).round()
    bins = xticks + 0.5
    xtickslabels = np.round(xs, 5)
    expL = []
    for c in range(1, 1 + maxCount):
        expC = SFS(c, mutNum, norm = False)
        expL.extend([c]*int(round(expC)))  ## get integer by round
    fig, ax = plt.subplots()
    ax.hist([expL, siteL], bins = bins, label = labelL, edgecolor = 'k', lw = 1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position(('axes', -0.01))
    plt.legend(loc = 'best')
    plt.xlabel('Site Frequency'); plt.ylabel('# of Sites')
    plt.xticks(bins[:-1], xtickslabels, rotation='vertical')
    plt.tight_layout()
    plt.savefig(outName, dpi = DPI, transparent = True)
    plt.close('all')

def Output(result,path,filename):
    with open("%s/%s.txt"%(path,filename),'w') as f:
        f.write("%s\t%s\t%s\n"%("cellpos_x","cellpos_y","cellmutation"))
        klist=list(result[0].keys())
        for k in klist:
            mutationL=traceLineage(result[0].mlineage,result[0][k].recentMut)
            mutationL.reverse()
            f.write("%s\t%s\t%s\n"%(k[0],k[1],','.join([str(x) for x in mutationL])))
    f.close()
 
# def compareNandS(b,d,m,maxSize):
#     args_N = []
#     parameter_N=[]
#     for i in np.arange(0,0.1,0.01):
#         initP = Population(initCells = [SpaceCell(b, d, m, '0', (0,0))], initT = 0, u = i, posR = 0, negR = 0, posS = 0, negS = 0, mlineage = None, initMutID = 0)
#         args_N.append((initP, maxSize))
#         parameter_N.append((b,d,m,i))
# 
#     args_S = []
#     parameter_S=[]
#     for i in np.arange(0,0.1,0.01):
#         for j in np.arange(0,0.1,0.01):
#             #s=np.random.exponential(0.1,1)[0]
#             s=0.1
#             initP = Population(initCells = [SpaceCell(b, d, m, '0', (0,0))], initT = 0, u = i, posR = j, negR = 0, posS = s, negS = 0, mlineage = None, initMutID = 0)
#             args_S.append((initP, maxSize))
#             parameter_S.append((b,d,m,i,j,s))
#     pools1 = multiprocessing.Pool(processes = 5)
#     pools2 = multiprocessing.Pool(processes = 5)
#     nL = pools1.starmap_async(rKMC, args_N)
#     sL = pools2.starmap_async(rKMC, args_S)
#     pools1.close(); pools1.join()
#     pools2.close(); pools2.join()
#     res_nL = [res for res in nL.get() if res]
#     res_sL = [res for res in sL.get() if res]
#     for i in range(len(res_nL)):
#         Output(res_nL[i],"./result","Netural_b-%s_d-%s_m-%s_u-%s"%parameter_N[i])
# 
#     for i in range(len(res_sL)):
#         Output(res_sL[i],"./result","Selection_b-%s_d-%s_m-%s_u-%s_pr%s_s-%s=0.1"%parameter_S[i])


def compareNandS(b,d,m,maxSize):
    args_N = []
    parameter_N=[]
    for i in np.arange(0,0.1,0.01):
        initP = Population(initCells = [SpaceCell(b, d, m, '0', (0,0))], 
                           initT = 0, 
                           u = i, 
                           posR = 0, negR = 0, 
                           posS = 0, negS = 0,
                           mlineage = None, 
                           initMutID = 0)
        args_N.append((initP, maxSize))
        parameter_N.append((b,d,m,i))

    pools1 = multiprocessing.Pool(processes = 10)
    nL = pools1.starmap_async(rKMC, args_N)
    pools1.close(); pools1.join()
    res_nL = [res for res in nL.get() if res]
    for i in range(len(res_nL)):
        Output(res_nL[i],"./res_simulation_mig_mut/","Netural_b-%s_d-%s_m-%s_u-%s"%parameter_N[i])

os.mkdir("./res_simulation_mig_mut")
if __name__ == '__main__':
  for j in np.arange(0, 0.2, 0.01):
    compareNandS(0.7, 0.1, j, 1e4)




    #startT = timer()
    #maxSize = 1e5
    #growthRate_compare(maxS = maxSize, reps = 3, proc = 6, outName = 'growthRate_compare_%.1e.jpg' % maxSize)
    #stopT = timer()
    #print('total runtime: %s' % (stopT - startT))
    
     # b, d, m = 0.793, 0.1, 0.1
     # reps, maxSize = 1, 2e4
     # startT = timer()
#     for i in range(reps):
#        initP = Population(initCells = [SpaceCell(b, d, m, '0', (0,0))], initT = 0, u = 0.1, posR = 0.1, negR = 0, posS = 0.1, negS = 0, mlineage = None, initMutID = 0)
#         ret = Gillespie(initP, maxSize)
#     gillT = timer()
#     print('totalTime (run Gillespie %s times to size %s): ' %(reps, maxSize), gillT - startT)
     # for i in range(reps):
         # initP = Population(initCells = [SpaceCell(b, d, m, '0', (0,0))], initT = 0, u = 0.1, posR = 0.1, negR = 0, posS = 0.1, negS = 0, mlineage = None, initMutID = 0)
         # ret = rKMC(initP, maxSize)
         # Output("./","ts.txt")     
#rKMCT = timer()  
    # print('totalTime (run rKMC %s times to size %s): ' %(reps, maxSize), rKMCT - gillT)
