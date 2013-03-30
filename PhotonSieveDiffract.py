#! usr/bin/env python
#-------------------------------------------------------------------------
#
#   File    :    PhotonSieveDiffract.py
#
#   Summary :    Class of functions to look at photon sieve diffraction
#                properties
#
#   Author  :    Robert Upton
#
#   Date    :    March 18, 2012, updated July 18, 2012, updated 08/17/2012
#
#-------------------------------------------------------------------------
import numpy as num
import pylab as py
import pdb
from matplotlib.patches import Circle
import scipy.special as sp

class sieveClass(object):
    '''This class contains the functions necessary to calculate the radii and
       locations of the holes in a photonsieve in a monolithic aperture (setHoleLocations)
       as well as their locations in a hexagonally segmented aperture (ZEMAXuserdefinedApertures,
        chooseSegment).
    
       The near-field and far field diffraction properties are calculated utilizing a 0th
       order Hankel transform that is centered on each photonsieve hole.  
            
        Ref:  Q. Cao and J. Jahns JOSA A. (2002)
        
        maxSDAp    : the maximum semi-diameter of the aperture
        fillFactor : the azimuthal duty cycle  
        rnVal, an, nHoles = setHoleLocations(maxSDap, masSDim, q, wavel, fillFactor) 
        
    '''

    def __init__(self,object):
        '''The initialization function that sets self as the object variable
        '''
        x = 1
        self.zPos = object
    
#-------------------------------------------------------------------------
    def setHoleLocations(self,maxSDap,q,wavel,fillFactor):
        '''The function that determines the locations of the diffraction holes in the photon sieve.  
        This function utilizes the phase matching conditions derived from Fresnel diffraction theory 
        to determine the locations of the Fresnel zones as a function of wavelength, and focal length.
           
           Ref:  Q. Cao and J. Jahns JOSA A. (2002)
            
           maxSDAp    : the maximum semi-diameter of the aperture
           q          : the focal length of the sieve
           wavel      : the wavelength of the radiation
           fillFactor : the azimuthal duty cycle  
           rnVal, an, nHoles = setHoleLocations(maxSDap, q, wavel, fillFactor) 
            
        '''
        maxorder   = num.floor((num.sqrt(maxSDap**2+q**2)-q)/wavel)     # the maximum radial order 
        
        rnVal      = num.zeros(maxorder)                                # the initialization of rnVal vector
        pi         = num.arccos(-1) 
        fillFactor = 1.0*fillFactor
        xp         = []
        yp         = []
        anIndx     = []

# determine radial locations of white zones for each order
        norder = 1
        while norder < maxorder:
            rnVal[norder] = num.sqrt((norder*wavel+q)**2-q**2)          # no binomial approx
            norder = norder + 1

        rnVal  = rnVal[num.nonzero(rnVal)]                              # empty rnVal elements = 0 
        nHoles = num.zeros(len(rnVal))                                  # initialize the nHoles vector    

# determine the sizes of the holes
# using the FRZ phase matching argument
        an      = num.zeros(len(rnVal))                                 # the vector containing the semi-diameter of holes
        wn      = num.zeros(len(rnVal))                                 # the widths of the clear zone
        circVal = num.zeros(len(rnVal))
        ll      = 0
        for jj in rnVal:
            wn[ll]      = wavel*q/2.0/jj                                # width of the positive zone
            an[ll]      = wn[ll]/2.0                                    # semi-diameter vector of holes
            circVal[ll] = 2*pi*jj
# circumference of zones
            nHoles[ll] = num.floor(circVal[ll]/fillFactor/2.0/an[ll])   # geometric model to determine number of holes
            for kk in num.arange(nHoles[ll]):                           
                dAngle = 2.0*pi*kk/nHoles[ll]                           # spacing of the angles    
                xp.append(rnVal[ll]*num.cos(dAngle))                    # positions of the holes    
                yp.append(rnVal[ll]*num.sin(dAngle))
            ll = ll + 1

# populate the self object
        self.an     = an
        self.rnVal  = rnVal
        self.nHoles = nHoles
        self.xp     = xp
        self.yp     = yp
        
        return self

#-------------------------------------------------------------------------
    def ZEMAXuserdefinedApertures(self,maxSDap,edgeD,Nrings,segPlotFlag):
        '''Generate the UDA files needed for ZEMAX to model a segmented photonsieve telescope.  
        This is only set-up for 2 rings right now (18 segments).  
        The hexagon segmented apertures are plotted if segPlotFlag == 1
        '''
        nHoles  = self.nHoles                   # nHoles vector elements are number of 
                                                # holes for each radial zone
        an      = self.an                       # an vector contains hole semi-diameter   
        xp      = self.xp                       # x locations of holes  
        yp      = self.yp                       # y locations of holes  
        pi      = num.arccos(-1)
        Nsegs   = 3*Nrings*(Nrings+1)           # number of hexgonal segments as function of rings  
        segDiam = 4.0*(maxSDap+edgeD)/5.0/num.sqrt(3)   # segment point-point diameter assuming 3rings
        segSep  = num.sqrt(3)/2.0*segDiam+edgeD
        xSeg    = num.zeros((18,1))             # segment x locations
        ySeg    = num.empty((18,1))             # segment y locations
        xPlotQ  = num.empty((0,))
        yPlotQ  = num.empty((0,))
        

        fidLocSeg = open('SegPosTable.txt','w')
# locate the positions of the segment centers        
        for segNum in num.arange(Nsegs):  
            if segNum <= 5:                                         # first ring
                xSeg[segNum] = segSep*num.sin(pi*segNum/3.0)
                ySeg[segNum] = segSep*num.cos(pi*segNum/3.0)
            elif segNum > 5 and segNum < 12:                        # second ring
                xSeg[segNum] = 2.0*segSep*num.sin(pi*segNum/3.0)
                ySeg[segNum] = 2.0*segSep*num.cos(pi*segNum/3.0)
            elif segNum > 11 & segNum < 19: 	                    # third ring			
                xSeg[segNum] = 2.0*segSep*num.sin(pi*segNum/3.0+pi/6.0)*num.cos(pi/6.0)
                ySeg[segNum] = 2.0*segSep*num.cos(pi*segNum/3.0+pi/6.0)*num.cos(pi/6.0)

# print the segment center locations to the file SegPosTable.txt
            strLoc = '%6.4f'%(xSeg[segNum]) + ' ' + '%6.4f'%(ySeg[segNum])   
            fidLocSeg.write('%s \n'%(strLoc))    
# choose plotting segments                  
            self.chooseSegment(segNum,segDiam,0,0,0,segPlotFlag) # determine location of holes in segment
            py.ion()
            py.figure(self.figSegNum,figsize=(11,11))
            py.plot(self.xPlotQ1 + xSeg[segNum],self.yPlotQ1 + ySeg[segNum],'k')    
            py.plot(self.xPlotQ2 + xSeg[segNum],self.yPlotQ2 + ySeg[segNum],'k')    
            py.plot(self.xPlotQ3 + xSeg[segNum],self.yPlotQ3 + ySeg[segNum],'k')    
            py.plot(self.xPlotQ4 + xSeg[segNum],self.yPlotQ4 + ySeg[segNum],'k') 
           
            fileStr="PrimSeg" + str(segNum+1) + ".UDA"
            fid    = open(fileStr,'w')                                          # generate a file handle for the given segment
            strVal = 'POL ' + '%6.4f '%(xSeg[segNum]) + '%6.4f '%(ySeg[segNum]) + \
                     ' ' + '%6.4f '%(segDiam/2.0) + '%i'%(6)
            fid.write('%s \n'%(strVal))                                         # write the segment 
        
            kk = 0
            for ii in num.arange(len(nHoles)):
                for jj in num.arange(nHoles[ii]):
                    dx = xp[kk] - xSeg[segNum]
                    dy = yp[kk] - ySeg[segNum]
                    self.chooseSegment(segNum,segDiam,dx,dy,an[ii],0)           # determine location of holes in segment
                    if self.segFlag == 1:                                       # if within segment
                        strVal = 'CIR ' + '%6.4f '%(xp[kk]) + '%6.4f '%(yp[kk]) + \
                                 '%6.4f'%(an[ii])
                        fid.write('%s \n'%(strVal))
                        py.plot(xp[kk],yp[kk])
                        ax   = py.gca()
                        xy   = [self.xp[kk],self.yp[kk]]
                        circ = Circle(xy, radius=self.an[ii])      # add circle plot 
                        ax.add_artist(circ)
                        circ.set_facecolor('none')     
                    kk = kk + 1                                    # update hole counter
            fid.close()                                            # close the .UDA file
       
        fidLocSeg.close()                                              # close the segment center location txt file        
        self.xSeg   = xSeg                                             # x location of segment
        self.ySeg   = ySeg                                             # y location of segment 
        self.xPlotQ = xPlotQ
        self.yPlotQ = yPlotQ

        return self

#-------------------------------------------------------------------------
    def chooseSegment(self,segNum,segDiam,dx,dy,holeSD,segPlotFlag):
        '''determine whether the hole is in the segment
           this code is based on a 4 quadrant model of an
           isomorphic hexagon.
           The hexagon is assumed to be circumscribed by
           a circle with radius SegDiam 
        '''
        segSD = segDiam/2.0
        if dx > 0:
            xHex  = dx + holeSD
        elif dx < 0:
            xHex = dx - holeSD
        else:
            xHex = dx
        
        if segPlotFlag:
            xPlotQ1 = num.empty((0,))
            yPlotQ1 = num.empty((0,))
            xPlotQ2 = num.empty((0,))
            yPlotQ2 = num.empty((0,))
            xPlotQ3 = num.empty((0,))
            yPlotQ3 = num.empty((0,))
            xPlotQ4 = num.empty((0,))
            yPlotQ4 = num.empty((0,))
        for ii in num.arange(4):
            Q = ii + 1
            if Q == 1:                              # quadrant 1
                yHex  = num.sqrt(3)*(segSD - xHex)
                xVals = [0, segSD/2.0, segSD]
                yVals = [0, num.sqrt(3)*segSD/2.0, yHex-holeSD]
                if segPlotFlag:
                    xPlotQ1 = num.append(xPlotQ1, [0, segSD/2.0, segSD])
                    yPlotQ1 = num.append(yPlotQ1, [num.sqrt(3)*segSD/2.0, num.sqrt(3)*segSD/2.0, 0])
                if xHex < xVals[1] and xHex > xVals[0] and dy > yVals[0] and dy <= yVals[1]:
                    self.segFlag = 1
                    break    
                elif xHex >= xVals[1] and xHex <= xVals[2] and dy >= yVals[0] and dy <= yVals[2]:
                    self.segFlag = 1
                    break    
                else:
                    self.segFlag = 0
            
            if Q == 2:                              # quadrant 2
                yHex  = num.sqrt(3)*(segSD + xHex)
                xVals = [-segSD/2.0, 0, -segSD]
                yVals = [0, num.sqrt(3)*segSD/2.0, yHex]
                if segPlotFlag:
                    xPlotQ2 = num.append(xPlotQ2, [0, -segSD/2.0,-segSD])
                    yPlotQ2 = num.append(yPlotQ2, [num.sqrt(3)*segSD/2.0, num.sqrt(3)*segSD/2.0, 0])
                if xHex < xVals[1] and xHex > xVals[0] and dy > yVals[0] and dy <= yVals[1]:
                    self.segFlag = 1
                    break    
                elif xHex <= xVals[0] and xHex >= xVals[2] and dy >= yVals[0] and dy <= yVals[2]:
                    self.segFlag = 1
                    break    
                else:
                    self.segFlag = 0
                        
            if Q == 3:                              # quadrant 3
                yHex  = -num.sqrt(3)*(segSD + xHex)
                xVals = [-segSD/2.0, 0, -segSD]
                yVals = [-num.sqrt(3)*segSD/2.0, 0, yHex]
                
                if segPlotFlag:
                    xPlotQ3 = num.append(xPlotQ3, [-segSD,-segSD/2.0, 0])
                    yPlotQ3 = num.append(yPlotQ3, [0, -num.sqrt(3)*segSD/2.0, -num.sqrt(3)*segSD/2.0])
                if xHex < xVals[1] and xHex > xVals[0] and dy > yVals[0] and dy <= yVals[1]:
                    self.segFlag = 1
                    break
                elif xHex <= xVals[0] and xHex >= xVals[2] and dy <= yVals[1] and dy >= yVals[2]:
                    self.segFlag = 1
                    break
                else:
                    self.segFlag = 0

            if Q == 4:                              # quadrant 4
                yHex  = -num.sqrt(3)*(segSD - xHex)
                xVals = [0, segSD/2.0, segSD]
                yVals = [-num.sqrt(3)*segSD/2.0, 0, yHex]
                
                if segPlotFlag:
                    xPlotQ4 = num.append(xPlotQ4, [0, segSD/2.0, segSD])
                    yPlotQ4 = num.append(yPlotQ4, [-num.sqrt(3)*segSD/2.0, -num.sqrt(3)*segSD/2.0, 0])
                if xHex < xVals[1] and xHex > xVals[0] and dy > yVals[0] and dy <= yVals[1]:
                    self.segFlag = 1
                    break
                elif xHex >= xVals[1] and xHex <= xVals[2] and dy <= yVals[1] and dy >= yVals[2]:
                    self.segFlag = 1
                    break
                else:
                    self.segFlag = 0

        # populate self with the plotting values of the segments   
        if segPlotFlag:
            self.xPlotQ1 = xPlotQ1
            self.xPlotQ2 = xPlotQ2
            self.xPlotQ3 = xPlotQ3
            self.xPlotQ4 = xPlotQ4
            self.yPlotQ1 = yPlotQ1
            self.yPlotQ2 = yPlotQ2
            self.yPlotQ3 = yPlotQ3 
            self.yPlotQ4 = yPlotQ4
            return self

#-------------------------------------------------------------------------
    def calcFarField(self,maxSDap,maxSDim,q,dq,FourierFlag,wavel,fillFactor,imSamp):
        '''Calculate the far-field component of the pinhole aperture
           Un(X,Y) = 2*Nf*An*exp(i*k*(Ln+r**2/2q))*J1(k*an*rho/q).
           The setHoleLocations function is called
           
           q : focal length  
            
        '''
# fundamental parameters
        pi = num.arccos(-1)
        k  = 2.0*pi/wavel
        
# call the aperture hole location function         
        self.setHoleLocations(maxSDap,q,wavel,fillFactor)
        
        xp     = self.xp
        yp     = self.yp
        nHoles = self.nHoles
        rnVal  = self.rnVal
        an     = self.an
        
# setup the aperture and image plane coordinates        
        xxi = num.matrix(num.linspace(-maxSDim,maxSDim,imSamp))
        yyi = num.transpose(xxi)
        unitVec = num.ones((imSamp,1))
        unitVec = num.matrix(unitVec)
        XX  = unitVec*xxi
        YY  = yyi*num.transpose(unitVec)

# calculate the complex field contributions        
        Un = num.zeros((imSamp,imSamp))
        kk = 0
        for ii in num.arange(len(nHoles)):
            Nf         = pi*an[ii]**2/q/wavel
            for jj in num.arange(nHoles[ii]):
                rho         = num.sqrt(num.power((XX-xp[kk]*(1-dq/q)),2)\
                            + num.power((YY-yp[kk]*(1-dq/q)),2))
                argBessel   = k*an[ii]*rho/q
                complexFac1 = Nf*num.exp(1j*k*num.power(rho,2)/2.0/q)
                complexFac2 = num.exp(1j*k*(xp[kk]**2+yp[kk]**2)*dq/2.0/q**2)
                complexFac  = complexFac1*complexFac2
                Jinc        = sp.jn(1,argBessel)/argBessel
                Un          = Un + num.multiply(complexFac,Jinc)
                kk          = kk + 1                 
        
        Utot = num.array(Un)
        Etot = num.conjugate(Utot)*Utot

# calculate the free space propagation to determine the focal region
# irradiance distribution
        if FourierFlag == 1:
            self.freeSpaceProp(Utot,xxi,dq,wavel,q,maxSDim,imSamp)
            Utotp = self.Utotp
        else:
            Utotp = Utot
            self.Utotp = Utotp
            
        Etot = num.conjugate(Utotp)*Utotp
        self.xxi  = xxi
        self.Etot = abs(Etot.real)
        self.Utot = self.Utotp
        
        return self

#-------------------------------------------------------------------------
    def calcNearField(self,maxSDap,maxSDim,q,dq,FourierFlag,wavel,fillFactor,imSamp):
        '''Calculate the far-field component of the pinhole aperture
            Un(X,Y) = 2*Nf*An*exp(i*k*(Ln+r**2/2q))*J1(k*an*rho/q).
            The setHoleLocations function is called
            
            q : focal length  
            
            '''
        # fundamental parameters
        pi = num.arccos(-1)
        k  = 2.0*pi/wavel
        
        # call the aperture hole location function         
        self.setHoleLocations(maxSDap,q,wavel,fillFactor)
        
        xp     = self.xp
        yp     = self.yp
        nHoles = self.nHoles
        rnVal  = self.rnVal
        an     = self.an
        
        # setup the aperture and image plane coordinates        
        xxi = num.matrix(num.linspace(-maxSDim,maxSDim,imSamp))
        yyi = num.transpose(xxi)
        unitVec = num.ones((imSamp,1))
        unitVec = num.matrix(unitVec)
        XX  = unitVec*xxi
        YY  = yyi*num.transpose(unitVec)
        
        # calculate the complex field contributions        
        Un = num.zeros((imSamp,imSamp))
        kk = 0
        for ii in num.arange(len(nHoles)):
            Nf         = pi*an[ii]**2/q/wavel
            for jj in num.arange(nHoles[ii]):
                rho         = num.sqrt(num.power((XX-xp[kk]*(1-dq/q)),2)\
                                       + num.power((YY-yp[kk]*(1-dq/q)),2))
                argBessel   = k*an[ii]*rho/q
                
                complexFac1 = 1j*Nf**2*num.exp(1j*k*num.power(rho,2)/2.0/q)
                complexFac2 = num.exp(1j*k*(xp[kk]**2+yp[kk]**2)*dq/2.0/q**2)
                complexFac  = (1-dq/q)*complexFac1*complexFac2
                Bessel      = sp.jn(1,argBessel)/argBessel-2*sp.jn(2,argBessel)/argBessel**2

                Un          = Un + num.multiply(complexFac,Bessel)
                kk          = kk + 1                 
        
        Utot = num.array(Un)
        Etot = num.conjugate(Utot)*Utot
        
        # calculate the free space propagation to determine the focal region
        # irradiance distribution
        if FourierFlag == 1:
            self.freeSpaceProp(Utot,xxi,dq,wavel,q,maxSDim,imSamp)
            UtotpNF = self.Utotp
        else:
            UtotpNF = Utot
            self.UtotpNF = UtotpNF
        
        EtotNF = num.conjugate(UtotpNF)*UtotpNF
        self.xxi  = xxi
        self.EtotNF = abs(EtotNF.real)
        self.UtotNF = self.UtotpNF
        
        return self

#-------------------------------------------------------------------------
    def calcFarFieldOnAxis(self,maxSDap,q,dq,Nplanes,wavel,fillFactor):
        '''calculate the on-axis field and irradiance
        '''

        # fundamental parameters
        pi = num.arccos(-1)
        k  = 2.0*pi/wavel
        
        # call the aperture hole location function
        self.setHoleLocations(maxSDap,q,wavel,fillFactor)
        
        xp     = self.xp
        yp     = self.yp
        nHoles = self.nHoles
        rnVal  = self.rnVal
        an     = self.an
        UtotZ  = []
        EtotZ  = []
        zPos   = []
        
        # calculate on-axis field
        for zz in num.arange(Nplanes):
            
            Un  = 0
            ddq = -dq/2.0 + zz*dq/Nplanes
            
            kk = 0
            for ii in num.arange(len(nHoles)):
                Nf         = pi*an[ii]**2/q/wavel
                for jj in num.arange(nHoles[ii]):
                    rho         = num.sqrt(xp[kk]**2+yp[kk]**2)*(1-ddq/q)
                    argBessel   = k*an[ii]*rho/q
                    complexFac1 = Nf*num.exp(1j*k*num.power(rho,2)/2.0/q)
                    complexFac2 = num.exp(1j*k*(xp[kk]**2+yp[kk]**2)*ddq/2.0/q**2)
                    complexFac  = complexFac1*complexFac2
                    Jinc        = sp.jn(1,argBessel)/argBessel
                    Un          = Un + num.multiply(complexFac,Jinc)
                    kk          = kk + 1
            
            zPos.append(ddq)
            Utot = num.array(Un)
            Etot = num.conjugate(Utot)*Utot
            UtotZ.append(Utot)
            EtotZ.append(abs(Etot.real))
        
        self.EtotOn = EtotZ
        self.UtotOn = UtotZ
        self.zPos   = zPos
        
        return self

#-------------------------------------------------------------------------
    def calcNearFieldOnAxis(self,maxSDap,q,dq,Nplanes,wavel,fillFactor):
        '''calculate the on-axis field and irradiance
        '''
        
        # fundamental parameters
        pi = num.arccos(-1)
        k  = 2.0*pi/wavel
        
        # call the aperture hole location function
        self.setHoleLocations(maxSDap,q,wavel,fillFactor)
        
        xp     = self.xp
        yp     = self.yp
        nHoles = self.nHoles
        rnVal  = self.rnVal
        an     = self.an
        UtotZ  = []
        EtotZ  = []
        zPos   = []
        
        # calculate on-axis field
        for zz in num.arange(Nplanes):
            
            Un  = 0
            ddq = -dq/2.0 + zz*dq/Nplanes
            kk = 0
            for ii in num.arange(len(nHoles)):
                Nf         = pi*an[ii]**2/q/wavel
                for jj in num.arange(nHoles[ii]):
                    rho         = num.sqrt(xp[kk]**2+yp[kk]**2)*(1-ddq/q)
                    argBessel   = k*an[ii]*rho/q
                    complexFac1 = 1j*Nf**2*num.exp(1j*k*num.power(rho,2)/2.0/q)
                    complexFac2 = num.exp(1j*k*(xp[kk]**2+yp[kk]**2)*ddq/2.0/q**2)
                    complexFac  = (1-ddq/q)*complexFac1*complexFac2
                    Bessel      = sp.jn(1,argBessel)/argBessel-2*sp.jn(2,argBessel)/argBessel**2
                    Un          = Un + num.multiply(complexFac,Bessel)
                    kk          = kk + 1
            
            zPos.append(ddq)
            Utot = num.array(Un)
            Etot = num.conjugate(Utot)*Utot
            UtotZ.append(Utot)
            EtotZ.append(abs(Etot.real))
        
        self.EtotNFOn = EtotZ
        self.UtotNFOn = UtotZ
        self.zNFPos   = zPos
        
        return self

#-------------------------------------------------------------------------
    def freeSpaceProp(self,Utot,xxi,dq,wavel,q,maxSDim,imSamp):
        '''free space Fresnel propagator
        '''
        
        dx   = 2.0*maxSDim/imSamp
        xi   = num.linspace(-0.5,0.5,imSamp)
        xi   = 1.0/dx*xi
        xi   = num.matrix(xi)
        eta  = num.transpose(xi)
        unitVec = num.ones((imSamp,1))
        unitVec = num.matrix(unitVec)
        Xi   = unitVec*xi
        Ei   = eta*num.transpose(unitVec)
        k    = 2*num.pi/wavel

# in the frequency domain
        UTOT  = num.fft.fft2(Utot)
        UTOT  = num.fft.fftshift(UTOT)
        Prop  = num.exp(-1j*k*(num.power(Xi,2)+num.power(Ei,2))*wavel**2*dq/2.0)
        UTOTp = num.multiply(UTOT,Prop)

# in the spatial domain
        Utotp = num.fft.ifft2(UTOTp)
        
        self.Utotp = Utotp
        
        return self

#-------------------------------------------------------------------------
    def ArrayCalcAlongZ(self,maxSDap,maxSDim,imSamp,wavel,fillFactor,path,fnamePNG,FourierFlag): 
        '''
        calculates the irradiance distribution along the optical axis of the instrument
        the 

        '''

        xxi        = num.linspace(-maxSDim,maxSDim,imSamp)
        zNum       = 0 
        
        for zz in self.zPos:
        
            zStep = zz
            ddq    = zStep
            print str(zStep)
            self.calcFarField(maxSDap,maxSDim,q,ddq,FourierFlag,wavel,fillFactor,imSamp)
            xp = self.xp
            yp = self.yp
            an = self.an
            nHoles = self.nHoles
            #self.calcNearField(maxSDap,maxSDim,q,ddq,FourierFlag,wavel,fillFactor,imSamp)
            IrradOn.append(self.Etot[imSamp/2.0,imSamp/2.0])
            zPlot.append(ddq)
                
# print hole locations and semi-diameters        
            if zz == 0:
                
                datHoles = num.zeros((len(self.xp),3))
                fid = open(fnameHoles,'w')
                cc = 0
                for ii in num.arange(len(nHoles)):
                    for jj in num.arange(nHoles[ii]):
                        datHoles[cc,0] = self.xp[cc]
                        datHoles[cc,1] = self.yp[cc]
                        datHoles[cc,2] = self.an[ii]
                        cc = cc + 1
                
                rr = 0    
                for rr in num.arange(len(datHoles)):
                    ss = 0
                    strVal = ''                 
                    for ss in num.arange(3):
                        strVal = strVal + '%12.6e '%(datHoles[rr,ss])
                    fid.write('%s\n'%(strVal))
                
                fid.close()
            
            # plot figures
            # aperture distribution    
            py.ion()
            py.figure(2,figsize=(12,4))
            py.clf()
            py.subplot(131)
            
            kk = 0
            for ii in num.arange(len(nHoles)):
                for jj in num.arange(nHoles[ii]):
                    py.plot(xp[kk], yp[kk])
                    ax   = py.gca()
                    xy   = [xp[kk],yp[kk]]
                    circ = Circle(xy, radius=an[ii])        # add circle plot 
                    ax.add_artist(circ)
                    circ.set_facecolor('none')
                    kk = kk + 1        
            
            py.xlabel(r'$\mathbf{X}(mm)$')
            py.ylabel(r'$\mathbf{Y}(mm)$')
            axis=py.gca()
            axis.axis('equal')
            py.ylim(-maxSDap,maxSDap)
            
            # profile
            #py.subplot(222)
            #py.plot(xxi,self.Etot[num.floor(imSamp/2.0),:]/num.max(self.Etot),'b')
            #py.plot(xxi,self.Etot[:,num.floor(imSamp/2.0)]/num.max(self.Etot),'k:')
            #py.title(r'$\mathbf{E}_{TOT}(\mathbf{x},0)$, 'r'$\mathbf{E}_{TOT}(0,\mathbf{y})$')
            #py.xlabel('position (mm)')
            #py.legend((r'$\mathbf{x}$',r'$\mathbf{y}$'),0)
            
            # total energy
            py.subplot(132)
            py.imshow(self.Etot,'bone')
            axis = py.gca()
            axis.axis('off')
            py.title(r'$\mathbf{E}_{TOT}$')
            
            # log of the total energy
            py.subplot(133)
            py.imshow(num.log10(self.Etot/num.max(self.Etot)),'bone')
            py.clim(-4,0)
            axis = py.gca()
            axis.axis('off')
            py.title(r'$log_{10}(\mathbf{E}_{TOT})$')
            
            figname = path + fnamePNG + str(zNum+1) + '.png'
            py.savefig(figname)
            zNum = zNum + 1
        
        # plot the on-axis irradiance as calculated by the full diffraction treatment
        #py.figure(1)
        #py.plot(zPlot,IrradOn/num.max(IrradOn),'kx')
        #py.legend(('Far-field','Near-field','TotalI','TotalII'),0)
        
        # produce a nice picture
        import ImageScienceToolBox as ISTB
        im = ISTB.ImageScienceClass()
        im.pol2cartMat(1000,1)
        RR   = im.RR
        mask = im.maskPol
        RR   = num.multiply(RR,mask)
        RR   = maxSDap*RR
        mask = RR
        FZmask = num.zeros((1000,1000))
        
        kk = 0
        for ii in self.rnVal:
            RR = (RR <= ii - self.an[kk]).choose(RR, 0)
            RR = (RR >= ii + self.an[kk]).choose(RR, 0)
            RR = (RR <> 0).choose(RR, 1)
            FZmask = FZmask + RR
            RR = mask
            kk = kk + 1
        
        py.figure(figsize = (10,5))
        py.subplot(121)
        py.imshow(FZmask, 'bone')
        axis = py.gca()
        axis.axis('off')
        py.title('Fresnel zones')
        
        py.subplot(122)
        kk = 0
        for ii in num.arange(len(nHoles)):
            for jj in num.arange(nHoles[ii]):
                py.plot(xp[kk], yp[kk])
                ax   = py.gca()
                xy   = [xp[kk],yp[kk]]
                circ = Circle(xy, radius=an[ii])        # add circle plot
                ax.add_artist(circ)
                circ.set_facecolor('none')
                kk = kk + 1
        
        axis = py.gca()
        axis.axis('equal')
        py.ylim(-maxSDap,maxSDap)
        py.title('Sieve holes')
        fname = path + 'NiceFZpic.png'
        py.savefig(fname)
        
        return self

#-------------------------------------------------------------------------
if __name__ == '__main__':
    
    maxSDap     = 20.0                                # semi-diameter of aperture
    q           = 20.0e03                              # focal length  
    edgeD       = 0.05                                 # separation between segments 
    Nrings      = 2                                    # number of hexagon rings 
    wavel       = 0.6328e-03                           # wavelength of radiation 
    Fratio      = q/2.0/maxSDap                        # Fratio for defocus parameter 
    maxSDim     = wavel*Fratio*30                      # image plane semi-diameter (mm) 
    FourierFlag = 0                                    # Fourier Flag = 0 no angular spectrum 
    fillFactor  = 2                                    # angular fill factor (form factor) 
    zPos        = [0]
    self        = sieveClass(zPos)                     # assign the class object to self
    self.segFlag = 1
    dq          = 0
    #dq          = 6.0*8.0*Fratio**2*wavel
    path        = '/users/rupton/documents/photonsieve/analysis/'
    imSamp      = 512
    
    pathHoles   = '/users/rupton/Documents/PHOTONSIEVE/ZEMAX/'
    fnameHoles  = pathHoles + 'PhHoles.txt'
    fnamePNG    = "photonSieveSBIR"
    self.figSegNum = 1                                 # figure number for plotting segments
    
    retVec = self.setHoleLocations(maxSDap,q,wavel,fillFactor)
    segPlotFlag = 1
    retVec = self.ZEMAXuserdefinedApertures(maxSDap,edgeD,Nrings,segPlotFlag)
    
    raw_input('%Press <return> to continue')