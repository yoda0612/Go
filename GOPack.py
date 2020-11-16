import numpy as np
import random


def cosCorner(z1,z2,z3):
    #[cosang] = cosCorner(z1,z2,z3) Get cos of angle at z1 in given triangle
    #Find cosine of angle at z1 in triangle with given centers
    l2=np.abs(z2-z1)
    l3=np.abs(z3-z1)
    l23=np.abs(z3-z2)
    denom=2*l2*l3
    cosang=(l2*l2+l3*l3-l23*l23)/denom
    cosang = max(min(cosang,1),-1)  # force to be between -1 and 1

    return cosang


class GOPack:
    def __init__(self):
        self.alpha=None
    def riffle(self,passNum):
        self.continueRiffle(passNum)

    def continueRiffle(self,passNum):
        cutval=0.01
        pass0=0
        maxVis = 2 * cutval
        while(pass0 < passNum and maxVis > cutval):
            self.layoutBdry()
            self.layoutCenters()
            self.setEffective()
            pass0+=1
            # print(pass0)
            break

    def setEffective(self):
        #interior adjustments
        for k in range(len(self.layoutVerts)):
            v = self.indx2v[k]
            targetArea = self.vAims[v]/2.0 # often negative for bdry vertices
            area = 0
            num = self.vNum[v]
            z = self.localcenters[v]
            flower = self.flowers[v]
            for j in range(num):
                jr = flower[j]
                jl = flower[j + 1]
                zr = self.localcenters[jr]
                zl = self.localcenters[jl]
                r = 0.5 * (np.abs(zr - z) + np.abs(zl - z) - np.abs(zr - zl))
                cC = cosCorner(self.localcenters[v], self.localcenters[jr], self.localcenters[jl])
                ang = np.arccos(cC)
                area = area + .5 * r * r * ang  # add area of sector

            # effective interior radius
            if (targetArea > 0.001):  #this one needs adjustment
                self.localradii[v] = (np.sqrt(area / targetArea)).real


        #bdry adjustments
        for k in range(self.bdryCount):
            w = self.bdryList[k]
            # targetArea related to aim
            targetArea = self.vAims[w] / 2.0 # often negative for bdry vertices
            angsum = 0
            area = 0
            num = self.vNum[w]
            z = self.localcenters[w]
            flower = self.flowers[w]
            for j in range(num):
                jr = flower[j]
                jl = flower[j + 1]
                if (self.bdryFlags[jr] >= 0 and self.bdryFlags[jl] >= 0): # neither is orphan
                    zr = self.localcenters[jr]
                    zl = self.localcenters[jl]
                    r = .5 * (np.abs(zr - z) + np.abs(zl - z) - np.abs(zr - zl))
                    cC = cosCorner(self.localcenters[w], self.localcenters[jr], self.localcenters[jl])
                    ang = np.arccos(cC)
                    angsum = angsum + ang
                    area = area + .5 * r * r * ang # add area of sector

            # effective radii
            if (targetArea > 0.001): # this one needs adjustment
                self.localradii[w] = (np.sqrt(area / targetArea)).real
                #the following is for (usually) bdry radii which would not be adjusted in traditional packing method, but are
                # adjusted in the GO approach.Change is averaged to moderate.
            elif (targetArea < -0.001): # e.g., bdry case in max packing
                self.localradii[w] = ((np.sqrt(2 * area / angsum)).real + self.localradii[w]) / 2.0



    def layoutCenters(self):
        self.updateVdata()

        layCount = len(self.layoutVerts)
        colCount = len(self.indx2v) - layCount
        #create sparse matrix A, m - by - m, where m = layCountt
        Aentries = np.zeros(self.tranIJcount,dtype=float)
        for k in range(self.tranIJcount):
            v = self.indx2v[self.tranI[k]]
            vrad = self.localradii[v]
            num = self.vNum[v]-1
            iR = self.inRadii[self.v2indx[v]]
            flower = self.flowers[v]
            j = self.tranJindx[k]
            if (j > -1):
                if (j == 0):
                    t1 = iR[num]
                else:
                    t1 = iR[j - 1]

                if (j > len(iR)):
                    continue
                t2 = iR[j]
                w = flower[j]
                Aentries[k] = ((t1 + t2) / (vrad + self.localradii[w])) / self.conduct[self.v2indx[v]]

            else: # edge to self if j=-1
                Aentries[k] = -1.0


        self.transMatrix = np.zeros((layCount, layCount),dtype=float)
        for i,xy in enumerate(zip(self.tranI,self.tranJ)):
            self.transMatrix[xy[0],xy[1]]=Aentries[i]

        # create matrix B, m x n where m = layCount, n = colCount
        Bentries = np.zeros(self.rhsIJcount, dtype=float)
        for k in range(self.rhsIJcount):
            v = self.indx2v[self.rhsI[k]]
            vrad = self.localradii[v]
            j = self.rhsJindx[k]
            num = self .vNum[v]-1
            iR = self.inRadii[self.v2indx[v]]
            flower = self.flowers[v]

            if (j == 0):
                t1 = iR[num]
            else:
                t1 = iR[j - 1]

            t2 = iR[j]
            w = flower[j]
            Bentries[k] = -1.0 * ((t1 + t2) / (vrad + self.localradii[w])) / self.conduct[self.v2indx[v]];

        #self.rhsMatrix = sparse(obj.rhsI, obj.rhsJ, Bentries, layCount, colCount);
        # self.rhsMatrix = np.zeros((self.rhsI, self.rhsJ), dtype=int)
        # Bentries_i = 0
        # for i in range(self.rhsI):
        #     for j in range(self.rhsJ):
        #         self.rhsMatrix[i, j] = Bentries[Bentries_i]
        #         Bentries_i += 1
        self.rhsMatrix = np.zeros((layCount, colCount), dtype=float)
        for i, xy in enumerate(zip(self.rhsI, self.rhsJ)):
            self.rhsMatrix[xy[0], xy[1]] = Bentries[i]

        # right hand side is B * vector of bdry centers
        zb = np.zeros(colCount, dtype=complex)
        for m in range(colCount):
            zb[m] = self.localcenters[self.indx2v[layCount + m]]
        self.rhs = np.dot(self.rhsMatrix,zb.reshape((-1,1))) # note the non-conjugate transpose

        #set centers
        Z=np.linalg.solve(self.transMatrix,self.rhs)
        for k in range(layCount):
            self.localcenters[self.indx2v[k]] = Z[k]
    def updateVdata(self):
        layCount = len(self.layoutVerts)
        if (layCount == 0):
            print('Error: "layoutVerts" was empty.\n')

        # update 'inRadii'
        self.inRadii = [None] * layCount
        for k in range(layCount):
            v = self.indx2v[k]
            vrad = self.localradii[v]
            flower = self.flowers[v]
            data = np.zeros(self.vNum[v],dtype=float)
            u = flower[0]
            urad = self.localradii[u]
            for j in range(self.vNum[v]):
                wrad = urad
                u = flower[j + 1]
                urad = self.localradii[u]
                data[j] = np.sqrt((vrad * urad * wrad) / (vrad + urad + wrad))
            self.inRadii[k] = data


        #update total conductances for 'layoutVerts'
        self.conduct = np.zeros(layCount,dtype=float)
        for k in range(layCount):
            v = self.indx2v[k]
            vrad = self.localradii[v]
            iR = self.inRadii[k]
            num = self.vNum[v]-1
            flower = self.flowers[v]
            w = flower[0]
            self.conduct[k] = (iR[num] + iR[0]) / (vrad + self.localradii[w])
            for j in range(1,num+1):
                t1 = iR[j - 1]
                t2 = iR[j]
                w = flower[j]
                self.conduct[k] = self.conduct[k] + (t1 + t2) / (vrad + self.localradii[w])

    def layoutBdry(self):
        self.setHoroCenters()

    def setHoroCenters(self):
        # initial for R: (sum of bdry radii) / pi
        R = 0.0
        minrad = 0.0
        r = np.zeros(self.bdryCount+1, dtype=float) # want closed list
        for j in range(self.bdryCount):
            r[j] = self.localradii[self.bdryList[j]]
            if (r[j] > minrad):
                minrad = r[j]
            R = R + r[j]
        r[self.bdryCount] = r[0] # close up
        R = R / np.pi
        if (R < 2.0 * minrad):
            R = 3.0 * minrad

        #Newton iteration to find R
        trys = 0
        keepon = 1
        while (keepon == 1 and trys < 100):
            trys = trys + 1
            fvalue = -2.0 * np.pi
            fprime = 0.0
            for j in range(self.bdryCount):
                Rrr = R - r[j] - r[j + 1]
                RRrr = R * Rrr
                ab = r[j] * r[j + 1]
                fvalue = fvalue + np.arccos((RRrr - ab) / (RRrr + ab))
                fprime = fprime - 1.0 * (R + Rrr) * np.sqrt(ab / RRrr) / (RRrr + ab)


            # is this working?
            newR = R - fvalue / fprime
            if (newR < R / 2.0):
                newR = R / 2

            if (newR > 2.0 * R):
                newR = 2.0 * R

            # cutoff(might be adjustable in future)
            if (np.abs(newR - R) < 0.00001):
                keepon = 0
            R = newR
        # end of while

        #scale all radii by 1 / R
        self.localradii = self.localradii/R;
        r = r/R;

        #set boundary centers, first being on y - axis.
        r2 = r[0]
        self.localcenters[self.bdryList[0]] = (1.0 - r2) *1j
        arg = np.pi/ 2.0;
        for k in range(1,self.bdryCount):
            r1 = r2
            r2 = r[k]
            RRrr = 1.0 - r1 - r2
            ab = r1 * r2
            delta = np.arccos((RRrr - ab) / (RRrr + ab))
            arg = arg + delta
            d = 1.0 - r2
            self.localcenters[self.bdryList[k]] = d * np.cos(arg) + 1j * d * np.sin(arg)



    def readpack(self,fname):
        data=[]
        with open(fname) as fp:
            for line in fp:
                if (len(line.strip())>0):
                    data.append(list(map(int, line.strip().split())))
        data=np.array(data)-1
        self.parse_triangles(data)
        self.indxMatrices()
        self.mode=1
        return self.nodeCount

    def parse_triangles(self,tList):
        self.edgeCount=0
        self.faceCount=tList.shape[0]
        bottom=tList.min()
        top=tList.max()

        #沒用
        #store initial count of faces containing each node
        _, clicks = np.unique(tList, return_counts=True)

        #make list of faces for each node
        nodefaces=[]
        for v in range(bottom,top+1):
            faces = np.argwhere(tList==v)[:,0]
            nodefaces.append(faces)


        #Start with first vertex of first face, build its flower
        tmpflower=[None]*(top+1)

        utilFlag = np.zeros((top+1), dtype=int) # -9, not yet encountered; -1, done (interior); -2, done (bdry); f, added by face f
        utilFlag.fill(-9)
        faceOrder = np.zeros(tList.shape[0], dtype=int) #1 ==> face touched (hence reoriented, if necessary)

        target = tList[0, 0]
        utilFlag[target] = 0
        faceOrder[0] = 1

        while (target != -1) :
            first_face = utilFlag[target]
            fvert = tList[first_face,:]
            ffindx = 0
            for i,fc in enumerate(nodefaces[target]):
                if (fc==first_face):
                    ffindx=i
                    break

            #make room, allowing for growth forward or back
            preflower = np.zeros(2 * nodefaces[target].shape[0] + 4 , dtype=int)

            #handle first face separately: put verts in middle of preflower space
            front = 0
            back = 0

            #put first two neighbors in middle and in proper order
            j,=np.where(fvert==target)
            v1 = fvert[(j+1)%3]
            v2 = fvert[(j+2)%3]

            back = nodefaces[target].shape[0]-1
            front = back + 1
            preflower[back] = v1;
            preflower[front] = v2;


            # indicate which face was used
            if (utilFlag[preflower[back]] == -9):
                utilFlag[preflower[back]] = first_face

            if (utilFlag[preflower[front]] == -9):
                utilFlag[preflower[front]] = first_face
            #print(utilFlag)

            nodefaces[target][ffindx] = -1 # this face has been used

            #now add petals forward(counterclockwise), then backward.
            #forward
            hit=1
            while (hit != 0 and (preflower[front] != preflower[back])):
                hit = 0
                v = preflower[front]

                for i,next_face in enumerate(nodefaces[target]):
                    if(next_face > -1 and hit == 0): # face not yet used
                        fvert = tList[next_face,:]
                        w = -1 # find edge(target, v) or (v, target)
                        j, = np.where(fvert == target)
                        va = fvert[(j + 1) % 3]
                        vb = fvert[(j + 2) % 3]

                        if (va == v):
                            w = vb
                            faceOrder[next_face] = 1
                        elif(vb == v): # must reverse this face
                            w = va
                            holdv = tList[next_face, 0]
                            tList[next_face, 0] = tList[next_face, 1]
                            tList[next_face, 1] = holdv
                            faceOrder[next_face] = 1
                        if (w > -1):
                            front = front + 1
                            preflower[front] = w;
                            if (utilFlag[w] == -9): # first encounter for w?
                                utilFlag[w] = next_face

                            hit = 1
                            nodefaces[target][i] = -1; # this face has been used

            #done with forward direction


            # flower open? must be bdry, so add petals backward
            if (preflower[front] != preflower[back]):
                hit = 1
                while (hit != 0 and preflower[front] != preflower[back]):
                    hit = 0
                    w = preflower[back]
                    for i,next_face in enumerate(nodefaces[target]):
                        if (hit == 0 and next_face > 0):
                            fvert = tList[next_face,:]
                            # find edge(target, w) or (w, target)
                            v = -1
                            j, = np.where(fvert == target)
                            wa = fvert[(j + 1) % 3]
                            wb = fvert[(j + 2) % 3]

                            if (wb == w):
                                v = wa
                                faceOrder[next_face] = 1
                            elif (wa == w): # must reverse this face
                                v = wb;
                                holdv = tList[next_face, 0]
                                tList[next_face, 0] = tList[next_face, 1]
                                tList[next_face, 1] = holdv
                                faceOrder[next_face] = 1;

                            if (v > -1):
                                back = back - 1
                                preflower[back] = v
                                if (utilFlag[v] == -9): # first encounter for v?
                                    utilFlag[v] = next_face
                                hit = 1
                                nodefaces[target][i] = -1;  # this face has been used
            #done with backward

            # done with target; fix up its tmpflower and mark as bdry / int
            tmpflower[target] = preflower[back:front+1]
            if (preflower[back] == preflower[front]):
                utilFlag[target] = -1 # interior
            else:
                utilFlag[target] = -2 # bdry

            # find the next vertex encountered but not done.
            target = -1;
            for v in range(top+1):
                if (utilFlag[v] > -1):
                    target = v;
                    break

            #break


        #Now 'tmpflower' should be complete. We need to check that
        newIndx = np.zeros((top+1),dtype=int)
        oldIndx = np.zeros((top+1),dtype=int)
        tick = -1
        nextv = 0
        while (nextv <= top):
            if (utilFlag[nextv] < 0): # % should have a flower for nextv
                tick = tick + 1;
                newIndx[nextv] = tick
                oldIndx[tick] = nextv

            nextv = nextv + 1
        oldIndx = oldIndx[0:tick+1]


        # #set nodeCount, allocate flowers, put reindexed verts in flowers
        self.nodeCount = (top+1)
        self.flowers = [None] * (top+1)
        self.vAims = np.zeros((top+1),dtype=float)
        self.vNum = np.zeros((top+1),dtype=int)
        bdryCount = 0
        nextv = 0
        while (nextv <= top):
            if (utilFlag[nextv] < 0):
                v = newIndx[nextv]
                n = len(tmpflower[nextv])
                self.vNum[v] = n - 1
                newflower = np.zeros(n, dtype=int)
                for i in range(n):
                    oldv = tmpflower[nextv][i];
                    newflower[i] = newIndx[oldv]
                self.flowers[v] = newflower;
                if (utilFlag[nextv] == -2): # bdry
                    bdryCount = bdryCount + 1;
                    self.gamma = v; # set gamma
                    self.vAims[v] = -1.0;
                else: #interior
                    self.vAims[v] = 2 * np.pi
            nextv = nextv + 1;

        #if alpha not already set, choose deep alpha
        if (self.alpha==None or self.alpha[0] == 0):
            seeds = np.zeros(self.nodeCount,dtype=int)
            tick = 0
            for j in range(self.nodeCount):
                if (utilFlag[j] == -2): #bdry
                    seeds[tick] = j;
                    tick = tick + 1;



            # find an alpha far from seed
            if (tick > 1):
                seeds = seeds[0:tick] # trim 'seeds'
                #print(seeds)

                alpha = self.FarVert(seeds)
                if (alpha < 0):
                    print('error in setting alpha\n')
                    self.alpha = 1
                else:
                    self.alpha = alpha

            else: # sphere? set 'alpha' to 1
                self.alpha = 1

        #organize combinatorics

        # organize combinatorics
        self.complex_count()
        if (bdryCount == 0):
            self.hes = 1 # sphere?
        else:
            self.hes = 0 # default to euclidean


        #finish with centers/radii. Note: calling routine needs to  call 'indxMatrices'
        facecount = self.faceCount

        self.radii = np.ones(self.nodeCount) * 0.5 # default to 1 / 2
        self.localradii = self.radii

        # store 'cents' if given, else default
        self.centers = np.zeros(self.nodeCount, dtype=complex)
        self.localcenters = self.centers



    def complex_count(self):
        #set 'tmpBdryFlags', 'alpha', and 'gamma'
        alpha = self.alpha
        flower = self.flowers[alpha]


        if (flower[0]!=flower[-1]):
            alpha = -1 # not interior, must reset

        gamma = self.gamma
        if (gamma == alpha):
            gamma = -1; # avoid setting to 'alpha'

        hasbdry = False
        tmpBdryFlags = np.zeros(self.nodeCount,dtype=int)
        for v in range(self.nodeCount):
            flower = self.flowers[v]
            if (flower[0] == flower[-1]): # interior?
                tmpBdryFlags[v] = 0
                if (alpha == -1):
                    alpha = v  # alpha defaults to first interior
            else:
                tmpBdryFlags[v] = 1
                if (gamma == -1):
                    gamma = v # gamma to first non - interior
                hasbdry = True

        if (alpha == -1):
            print('Error: complex has no interior vertex')
            nodecount = -1
            return nodecount

        self.alpha = alpha

        #sort and mark vertices
        self.bdryList = []
        self.layoutVerts = []
        self.orphanVerts = []
        self.orphanEdges = []

        #spherical case
        if (not hasbdry):
            self.hes = 1
            # cclw faux bdry {a, b, c}, some face far from alpha
            a = self.FarVert([self.alpha])
            flower = self.flowers[a]
            b = flower[1]
            c = flower[0]
            self.bdryList = [a, b, c, a]
            self.gamma = a
            self.intVerts = np.zeros(self.nodeCount - 3,dtype=int)
            tick = 0
            for v in range(self.nodeCount):
                if (v!=a and v!=b and v!=c):
                    self.intVerts[tick] = v
                    tick = tick + 1;


            self.orphanVerts = [];
        else:
            #non-spherical
            status = np.zeros(self.nodeCount,dtype=int)
            hitlist=[alpha]
            status[alpha] = -1
            self.intVerts = []
            self.intVerts.append(alpha)

            while (len(hitlist)!=0):
                # pick off the next vert to process
                v = hitlist[0]
                hitlist = hitlist[1:]
                flower = self.flowers[v]
                for k in range(len(flower)):
                    w = flower[k]
                    if (tmpBdryFlags[w] == 0 and status[w] == 0): # interior, not hit before
                        status[w] = -1
                        hitlist.append(w)
                        self.intVerts.append(w)

                status[v] = 1 # finished with v

            # create 'tmpBdryV' list
            # Set status(w) = -1 if w neighbors interior vert
            lolong = len(self.intVerts)
            tmpBdryV = []
            for j in range(lolong):
                flower = self.flowers[self.intVerts[j]]
                for k in range(len(flower)):
                    w = flower[k]
                    if (status[w] == 0): # w not yet encountered
                        status[w] = -1
                        tmpBdryV.append(w)

            # Organize 'bdryList' so it is cclw order
            # try to start with 'gamma', else set 'gamma'
            self.bdryList = np.zeros(len(tmpBdryV) + 1,dtype=int)
            firstbdry = self.gamma;
            if (firstbdry < 1 or firstbdry > self.nodeCount or status[firstbdry]!=-1): # choose new start
                firstbdry = tmpBdryV[0] # first bdry encountered
                self.gamma = firstbdry

            # fill 'bdryList'; next nghb is first petal with contact
            tick = 0
            self.bdryList[tick] = firstbdry

            #start with 'firstbdry'
            flower = self.flowers[firstbdry]
            for j in range(len(flower)):
                nextb = flower[j]
                if (status[nextb] == -1):
                    tick = tick + 1
                    self.bdryList[tick] = nextb
                    break

            # cycle through from 'nextb'
            while (nextb!=firstbdry and tick < 2 * len(tmpBdryV)):
                flower = self.flowers[nextb]
                for j in range(len(flower)):
                    nextb = flower[j]
                    if (status[nextb] == -1): # first downstream neighboring interior
                        tick = tick + 1
                        self.bdryList[tick] = nextb
                        break

            if (nextb!=firstbdry or self.bdryList[-1]!=self.bdryList[0]):
                print('Error forming bdryList')

            if (tick == 2 * len(tmpBdryV)):
                print('Error in bdryList, too long')

            # are there orphans?
            for i in range(self.nodeCount):
                if (status[i] == 0):
                    self.orphanVerts.append(i)

        #set the counts
        nodecount = self.nodeCount
        self.intCount = len(self.intVerts)
        self.bdryCount = len(self.bdryList) - 1
        self.orphanCount = len(self.orphanVerts)

        self.edgeCount = 0
        self.faceCount = 0
        self.vNum = np.zeros(self.nodeCount,dtype=int)
        self.bdryFlags = np.zeros(self.nodeCount,dtype=int)
        totNum = 0

        for v in range(self.nodeCount):
            flower = self.flowers[v]
            num = len(flower) - 1
            self.vNum[v] = num
            totNum = totNum + self.vNum[v]
            for k in range(num):
                w = flower[k]
                if (w > v):
                    self.edgeCount = self.edgeCount + 1
            if (flower[0]!=flower[-1]):
                w = flower[-1]
                if (w > v):
                    self.edgeCount = self.edgeCount + 1;
                self.bdryFlags[w] = 1
        self.faceCount = totNum / 3

        #set up orphanEdges
        if (len(self.orphanVerts)!=0):
            self.orphanEdges = np.zeros(self.orphanCount, 2)

            # keep track with ctlg = -2 for interior, -1 for bdryList,
            ctlg = np.zeros(self.nodeCount,dtype=int)
            for j in range(len(self.intVerts)):
                v = self.intVerts[j]
                ctlg[v] = -2
            for j in range(len(self.bdryList) - 1):
                w = self.bdryList[j]
                ctlg[w] = -1

            #find verts next to just two bdry verts, catalog them and store their bdry edges.
            tick = 0
            for j in range(self.orphanCount):
                v = self.orphanVerts[j]
                flower = self.flowers[v]
                for k in range(self.vNum[v]):
                    m = flower[k]
                    n = flower[k + 1]
                    if (ctlg[m] == -1 and ctlg[n] == -1):  # successive bdryList nghbs
                        self.orphanEdges[j,:]=[n, m]
                        ctlg[v]= j
                        tick = tick + 1;
                        break


            #rest eventually next to catalogued ones and inherit their bdry edges
            hit = 1
            while (hit > 0):
                hit = 0
                for j in range(self.orphanCount):
                    v = self.orphanVerts[j]
                    flower = self.flowers[v]

                    if (ctlg[v] <= 0):
                        for k in range(len(flower)):
                            w = flower[k]
                            if (ctlg(w) > 0): # neighbor is catalogued orphan
                                hit = hit + 1
                                tick = tick + 1
                                ctlg[v] = ctlg[w]
                                self.orphanEdges[j,:]=self.orphanEdges[ctlg[w],:]

            if (tick < self.orphanCount):
                print('Error: seem to have missed some orphans');

        # default 'layoutVerts', 'rimVerts'
        self.layoutVerts = self.intVerts;
        self.rimVerts = self.bdryList;
        self.v2indx = [];
        self.indx2v = [];


    def FarVert(self,seeds):
        nodecount = self.nodeCount
        marks = np.zeros((nodecount, 1),dtype=int)

        #juggle two lists, 'curr' and 'nextlist'
        nextlist = seeds
        gennum = 1
        farvert = seeds[0]
        nothits = nodecount
        while (nextlist.shape[0]!=0 and (gennum < 10 or nothits < 100)):
            # debug fprintf('gennum=%d\n', gennum);
            getnonz = 0
            for j in range(nextlist.shape[0]):
                if (nextlist[j] > 0):
                    getnonz = getnonz + 1


            curr = nextlist[0:getnonz+1]
            nextlist = np.zeros(getnonz * 5,dtype=int)
            nextend = 0
            for j in range(curr.shape[0]):
                k = curr[j]
                if(marks[k] < 1):
                    nothits = nothits - 1

                marks[k] = gennum
                farvert = k
                #print(k)
                flower = self.flowers[k]
                for m in range(len(flower)):
                    p = flower[m]
                    if (marks[p] == 0):
                        nextlist[nextend ] = p
                        nextend = nextend + 1

            nextlist = nextlist[0:nextend] # trim zeros
            #print("nextlist",nextlist)
            gennum = gennum + 1

        #if we visited most vertices, current 'farvert' should be deep
        if (nothits <= 100):
            return farvert

        # reaching here, should have searched to depth 10 with lots
        farvert=0
        while (farvert==0):
            k=random.randint(0,nodecount-1)
            if (marks[k]==0):
                farvert=k
                return farvert

    def indxMatrices(self):
        self.layoutVerts = self.intVerts
        self.rimVerts = self.bdryList

        lolong = len(self.layoutVerts)
        #build 'indx2v' and 'v2indx'; 'layoutVerts' first, then 'rimVerts'
        self.v2indx = np.zeros(self.nodeCount, dtype=int)
        self.indx2v = []
        for j in range(lolong):
            self.indx2v.append(self.layoutVerts[j])
            self.v2indx[self.layoutVerts[j]] = j

        # and finish with 'rimVerts'
        for j in range(len(self.rimVerts) - 1):
            w = self.rimVerts[j]
            self.indx2v.append(w)
            self.v2indx[w] = len(self.indx2v)

        #set up tranI / J / Jindx, rhsI / J / Jindx; cut down to right size later
        ijCount = 0
        for k in range(lolong):
            ijCount = ijCount + self.vNum[self.indx2v[k]] + 1

        tranI = np.zeros(ijCount, dtype=int)
        tranJ = np.zeros(ijCount, dtype=int)
        tranJindx = np.zeros(ijCount, dtype=int)
        rhsI = np.zeros(ijCount, dtype=int)
        rhsJ = np.zeros(ijCount, dtype=int)
        rhsJindx = np.zeros(ijCount, dtype=int)

        kj = 1  # count the tran entries
        kw = 1  # count the rhs entries
        for k in range(lolong):
            v = self.indx2v[k]

            # diagonal entry first
            tranI[kj] = k
            tranJ[kj] = k
            tranJindx[kj] = -1 # indicates edge to self
            kj = kj + 1

            # now petal entries
            flower = self.flowers[v]
            num = len(flower) - 1
            for j in range(num):
                w = flower[j]
                if (self.v2indx[w] <= lolong): #interior petal
                    tranI[kj] = k
                    tranJ[kj] = self.v2indx[w]
                    tranJindx[kj] = j
                    kj = kj + 1
                else: # bdry petal
                    rhsI[kw] = k
                    rhsJ[kw] = self.v2indx[w] - lolong # offset by layCount for later use
                    rhsJindx[kw] = j
                    kw = kw + 1


        # store this info at the proper sizes
        self.tranIJcount = kj - 1
        self.rhsIJcount = kw - 1
        self.tranI = tranI[1:kj]
        self.tranJ = tranJ[1:kj]
        self.tranJindx = tranJindx[1:kj]
        self.rhsI = rhsI[1:kw]
        self.rhsJ = rhsJ[1:kw]-1
        self.rhsJindx = rhsJindx[1:kw]

    def __str__(self):
        #print("alpha",self.alpha)
        # print("gamma",self.gamma)
        # print("nodeCount",self.nodeCount)
        # print("edgeCount",self.edgeCount)
        # print("faceCount",self.faceCount)
        #print("vNum",self.vNum)
        # print("flower",self.flowers)
        # print("bdrtFlags",self.bdryFlags)
        # print("vAims",self.vAims)
        # print("intVerts",self.intVerts)
        # print("bdryList",self.bdryList)
        # print("layoutVerts",self.layoutVerts)
        # print("rimVerts",self.rimVerts)
        # print("tranI",len(self.tranI),self.tranI)
        # print("tranJ",len(self.tranI),self.tranJ)
        # print("tranJindx",len(self.tranJindx),self.tranJindx)
        # print("rhsI", len(self.rhsI), self.rhsI)
        # print("rhsJ", len(self.rhsJ), self.rhsJ)
        # print("rhsJindx", len(self.rhsJindx), self.rhsJindx)
        #print("localradii",self.localradii)

        return ""