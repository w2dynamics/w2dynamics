import numpy as np
import scipy.integrate
import MAXENT

def kkt(aws, ws, integrator=np.trapz):
    """ Does a Kramers-Kronig transform for the spectral function:
                Re G(w) = p.v. integral  A(w')/(w - w') dw'
    """
    realgw = np.empty_like(aws)
    for iw, w in enumerate(ws):
      denom = w - ws
      denom[iw] = 1  # to avoid division by zero
      integrand = (aws - aws[iw])/denom
      realgw[iw] = integrator(integrand, ws)
    return realgw
 
class MaxEntBase:
   """ 
   Base class of the python wrapper to Benedikts maxent code
   this class defines all parameters accesible to the user and
   sets some default values for the Monte Carlo procedure
      
   further simple methods are included to allocate arrays 
   (self.allocateArrays) and compute the omega grid (self.computeGrid)
      
   the f2py call is wrapped by the method self.computeSpec
   a more detailed documentation can be found in the fortran program
   itself and in Benedikts project report "Maximum Entropy Method for 
   Quantum Monte Carlo Simulations"
  
   Author: Benedikt Hartl
   Author: Patrik Gunacker
   """
      
   def __init__(self,kernelMode,beta,Ncorr,NGrid1,NGrid2,NGrid3,\
                     wmin,w1,w2,wmax,NSmooth,smoothCnt,smoothingMethod,logprefix):
      #possible kernels for the analytic continuation:
      #0-> gtau->A(w) 
      #1-> giw->A(w)
      #10 -> chi(tau)->A(w)=chi(w)/w
      #11 -> chi(tau)->A(w)=chi(w)/(1-exp(...))
      self.kernelMode=kernelMode
      
      #inverse temperature
      self.beta = beta

      #number of data points of the Greens/Chi input function
      self.Ncorr = Ncorr
      
      #number of grid points in region 1,2,3
      self.NGrid1 = NGrid1
      self.NGrid2 = NGrid2
      self.NGrid3 = NGrid3
      self.NGrid = NGrid1 + NGrid2 + NGrid3
      
      #frequency range in region 1,2,3
      self.wmin=wmin
      self.w1=w1
      self.w2=w2
      self.wmax=wmax
      
      #smoothing parameters: NSmooth -> number of smoothing steps
      #smoothCnt -> number of neighbors to smooth
      self.NSmooth = NSmooth
      self.smoothCnt = smoothCnt
      
      #smoothing method
      #0: locally average spectral function: 
      #   smoothSpectrum(omega_i)=sum(spectralFunction(k:l))/(l-k+1), k=omega_i-smoothCnt, l=omega_i+smoothCnt, 
      #1: smoothing by locally integrate spectral function: 
      #   smoothSpectrum(omega_i)=SUM(domega(k:l)*spectralFunction(k:l))/SUM(domgea(k:l), important for non-equidistant grid-points 
      #2: smoothing by locally integrate chi''(omega) propTo spectral function
      #   important for non-equidistand grid points for kernel 10 and 11
      self.smoothingMethod = smoothingMethod
      
      #prefix for log files
      self.logprefix = logprefix
      
      #setting some sane default parameters
      #default entropy weight
      self.alpha = 1000

      #monte carlo steps in the variation of maxent
      self.NMonteCarlo = 6000
      
      #annealing steps in the variation of maxent
      self.NAnnealing = 100
         
      #default seed
      self.seed = -1431955233
      
      #checking for range sane input and allocating arrays
      self.checkRange()
      self.allocateArrays()
            
   
   def checkRange(self):
      #checking frequency range input for order
      if(self.wmin > self.w1):
         raise ValueError("wmin > w1")
      if(self.w1 > self.w2):
         raise ValueError("w1 > w2")
      if(self.w2 > self.wmax):
         raise ValueError("w2 > wmax")         
   
   def computeGrid(self):
      #creating w-grid for the three regions
      self.grid=np.zeros((self.NGrid),dtype=np.float)
      self.grid[:self.NGrid1]=np.linspace(self.wmin,self.w1,num=self.NGrid1,endpoint=False)
      self.grid[self.NGrid1:self.NGrid1+self.NGrid2]=np.linspace(self.w1,self.w2,num=self.NGrid2,endpoint=True)
      self.grid[-self.NGrid3:]=np.linspace(self.wmax,self.w2,num=self.NGrid3,endpoint=False)[::-1]

   def allocateArrays(self):        
      #arrays in input data points
      self.corrGrid = np.zeros(self.Ncorr,dtype=float,order='F')
      self.corrCov = np.zeros(self.Ncorr,dtype=float,order='F')
      self.corrReal = np.zeros(self.Ncorr,dtype=float,order='F')
      self.corrImag = np.zeros(self.Ncorr,dtype=float,order='F')
      self.corr = np.zeros(self.Ncorr,dtype=complex,order='F')
      
      #arrays in output points
      self.grid = np.zeros(self.NGrid,dtype=float,order='F')
      self.model = np.zeros(self.NGrid,dtype=float,order='F')
      self.spec = np.zeros(self.NGrid,dtype=float,order='F')
      self.gws = np.zeros(self.NGrid,dtype=np.complex,order='F')

   def computeSpec(self):
      #f2py call
      MAXENT.mmaximumentropy.w2maxent (self.spec , self.model, self.grid, self.corr, self.corrCov, self.corrGrid,\
               self.alpha, self.kernelMode, self.beta, self.NMonteCarlo, self.NAnnealing, self.seed, self.NSmooth , self.smoothCnt,\
               self.smoothingMethod,self.logprefix)

class Green2Aw(MaxEntBase):
   """
   this derived class is used for the continuation of a greens function in tau (kernel mode 0) and iw (kernel mode 1)
   some default values for the output window and the smoothing method are set
   the model is assumed to be flat and normalized
   the ouput is adapted for the type of continuation
   """
   def __init__(self,kernelMode,beta,\
                     Ncorr,NGrid1=100,NGrid2=101,NGrid3=100,\
                     wmin=-30.,w1=-10.,w2=10.,wmax=30.,\
                     NSmooth=3,smoothCnt=3,smoothingMethod=0,logprefix="out"):
      
      if(kernelMode != 0 and kernelMode != 1):
         raise ValueError("Please choose valid kernel 0 or 1")
      
      MaxEntBase.__init__(self,kernelMode,beta,Ncorr,NGrid1,NGrid2,NGrid3,\
                          wmin,w1,w2,wmax,NSmooth,smoothCnt,smoothingMethod,logprefix)
   
   def symmetrize(self):
      
      self.corrReal[:self.corrReal.shape[0]//2] = 0.5*(self.corrReal[:self.corrReal.shape[0]//2] 
                                                  + self.corrReal[self.corrReal.shape[0]//2:]) 
      self.corrReal[self.corrReal.shape[0]//2:] = self.corrReal[:self.corrReal.shape[0]//2]
     
      self.corrCov[:self.corrCov.shape[0]//2] = 0.5*(self.corrCov[:self.corrCov.shape[0]//2] 
                                                  + self.corrCov[self.corrCov.shape[0]//2:]) 
      self.corrCov[self.corrCov.shape[0]//2:] = self.corrCov[:self.corrCov.shape[0]//2]
  
   def computeModel(self):
      self.model[:] = 1.
      MAXENT.mmaximumentropy.normalizetogrid(self.model,self.grid,self.kernelMode)
   
   def computeGws(self):
      """Takes A(w) and computes the (complex) G(w) using a KKT."""
      self.gws.real = kkt(self.spec, self.grid)
      self.gws.imag = -np.pi*self.spec
      
   def printSpec(self,fname):
      aw = file(fname, "w")  
      print >> aw,  "# w      A(w)      model(w)"
      for i in xrange(self.NGrid):
         print >> aw, self.grid[i], self.spec[i], self.model[i]
         
   def printGws(self, fname):
      gws = file(fname, "w")  
      print >> gws,  "# w      Re[G(w)]      Im[G(w)]      model(w)"
      for i in xrange(self.NGrid):
         print >> gws, self.grid[i], self.gws[i].real, self.gws[i].imag, self.model[i]
         
       
class Chi2Aw(MaxEntBase):
   """
   this derived class is used for the continuation of bosonic chi function in tau (kernel mode 10 or 11)
   kernel mode 10 extends the spectral function A(w) by a factor 1/w
   kernel mode 11 extends the spectral function A(w) by a factor 1/exp(...)
   chi(w) is calculated from A(w) during the output including a factor pi
   some default values for the output window and the smoothing method are set
   the model is assumed as a half-gaussian step
   """
   def __init__(self,kernelMode,beta,\
                     Ncorr,NGrid1=120,NGrid2=60,NGrid3=30,\
                     wmin=0.,w1=10.,w2=20.,wmax=40.,\
                     NSmooth=3,smoothCnt=3,smoothingMethod=2,logprefix="out"):
      
      if(kernelMode != 10 and kernelMode != 11):
         raise ValueError("Please choose valid kernel 10 or 11")
      
      MaxEntBase.__init__(self,kernelMode,beta,Ncorr,NGrid1,NGrid2,NGrid3,\
                          wmin,w1,w2,wmax,NSmooth,smoothCnt,smoothingMethod,logprefix)
      
   """
   some futher possible default values
   NGrid1=120,NGrid2=60,NGrid3=15       # kernel 10
   wmin=0,w1=10D0,w2=20D0,wmax=30D0     # kernel 10
   NGrid1=120,NGrid2=60,NGrid3=30       # kernel 11
   wmin=0 , w1=10D0,w2=20D0,wmax=40D0   # kernel 11
   """
      
   def computeModel(self):
      #FIXME: at some point we replace this by a more physically motivated model
      sigmaModel=1.
      heightModel=0.4
      for i in xrange(self.NGrid):
         self.model[i]=np.exp(-self.grid[i]**2/(2./sigmaModel**2))
         #model[i]=1./(1.+np.exp((grid[i]-sigmaModel)*5))

      self.model[:]=self.model[:]/self.model[1]
      self.model[:]=self.model[:]*heightModel  
   
   def computeGws(self):
      """Takes A(w) and computes the (complex) Chi(w) using a KKT."""
      
      kernelModeFactor = np.ones(self.NGrid,dtype=float,order='F')
      
      for i in xrange(self.NGrid):
         if self.kernelMode == 10: 
             kernelModeFactor[i] = self.grid[i]*np.pi

         elif self.kernelMode == 11:
             kernelModeFactor[i] = (1 - np.exp(-self.beta*self.grid[i]))*np.pi
      
      #FIXME: this is not yet tested for X(w) real and imaginary part
      self.gws.real = kkt(kernelModeFactor*self.spec, self.grid)
      self.gws.imag = kernelModeFactor*self.spec
   
   def printSpec(self,fname):
      aw = file(fname, "w")
      if self.kernelMode == 10:
         print >> aw,  "# w   Chi(w)   Chi(w)/(w*pi)    model(w)"
      elif self.kernelMode == 11:
         print >> aw,  "# w   Chi(w)   Chi(w)/(exp(...)*pi)    model(w)"
         
      for i in xrange(self.NGrid):
         kernelModeFactor = 1
         
         #bosonic kernel -> chi/w from max ent
         if self.kernelMode == 10: 
            kernelModeFactor=self.grid[i]*np.pi
            #if (i == 1):
            #   kernelModeFactor = 1e-12
         
         #bosonic kernel -> chi/bose from max ent
         elif self.kernelMode == 11:
            kernelModeFactor = (1 - np.exp(-self.beta*self.grid[i]))*np.pi
            #    if (i .EQ. 1):
            #       kernelModeFactor = 1e-12
         
         print >> aw, self.grid[i], self.spec[i]*kernelModeFactor, self.spec[i], self.model[i]

   def printGws(self, fname):
      gws = file(fname, "w")
      print >> gws,  "# w      Re[Chi(w)]      Im[Chi(w)]      model(w)"
      for i in xrange(self.NGrid):
         print >> gws, self.grid[i], self.gws[i].real, self.gws[i].imag, self.model[i]

