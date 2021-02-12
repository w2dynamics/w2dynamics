from __future__ import absolute_import, print_function, unicode_literals
import numpy as np
try:
    from scipy.special import factorial
except ImportError:
    from scipy.misc import factorial
import re
import sys
import w2dyn.dmft.dynamicalU as dynamicalU

def udensity_values(nbands=1, u=0., v=0., j=0.):
    """Creates the U matrix for the density-density interaction"""
    udens = np.zeros((nbands, 2, nbands, 2))
    if nbands > 1:
        if v: udens[...] = v
        if j: udens -= np.reshape((j, 0, 0, j), (1, 2, 1, 2))
    # the norbitals > 1 statements also set diagonal elements,  so we
    # unconditionally need to set them here
    band = np.arange(nbands)
    udens[band, :, band, :] = np.array(((0, u), (u, 0)))
    return udens

def screenedudensity_values(nbands=1, Uw_Mat=0, u=0., v=0., j=0.):
    """Creates the screened U matrix for the density-density interaction"""
    udens = np.zeros((nbands, 2, nbands, 2))
    zero_umatrix = np.zeros((nbands*2,nbands*2,nbands*2,nbands*2))
    _Retarded2Shifts = dynamicalU.Replace_Retarded_Interaction_by_a_Shift_of_instantaneous_Potentials(0, zero_umatrix, Uw_Mat, 1)
    
    for i in range(0,nbands):
      udens[i,1,i,0] = u + _Retarded2Shifts.U_shift[i]
      udens[i,0,i,1] = u + _Retarded2Shifts.U_shift[i]

    n = 0
    for i in range(0,nbands):
      for l in range(i+1,nbands):
        if (_Retarded2Shifts.Vw is True):
          v_s = v + _Retarded2Shifts.V_shift[n]
        else:
          v_s = v
        udens[i,:,l,:] = v_s
        udens[l,:,i,:] = v_s
        n+=1

    for i in range(0,nbands):
      for l in range(i+1,nbands):
        if (_Retarded2Shifts.Jw is True):
          j_s = j + _Retarded2Shifts.J_shift[n]
        else:
          j_s = j
        udens[i,1,l,1] -= j_s
        udens[i,0,l,0] -= j_s
        udens[l,1,i,1] -= j_s
        udens[l,0,i,0] -= j_s
        n+=1
    return udens

def udensity_from_ufull(ufull):
    r"""Extracts  the density-density part from a full U-matrix

    .. math::   U_{i\sigma,j\tau} = U_{ijij} - U_{ijji} \delta_{\sigma\tau}

    Example:
    --------
    >>> ufull = ufull_kanamori(2, 10.0, 1.0, 0.1)
    >>> udens = udensity_values(2, 10.0, 1.0, 0.1)
    >>> np.allclose(udens, udensity_from_ufull(ufull))
    True
    """
    nbands = ufull.shape[0]
    spdiag = 0, 1
    udens = np.zeros((nbands, 2, nbands, 2), ufull.dtype)
    udens[...] = ufull.diagonal(0, 0, 2).diagonal(0, 0, 1)[:,None,:,None]
    udens[:,spdiag,:,spdiag] -= ufull.diagonal(0, 1, 2).diagonal(0, 0, 1)
    return udens

def udensity_from_spindep_ufull(ufull):
    r"""Extracts  the density-density part from a full U-matrix
    .. math::   U_{i\sigma,j\tau} = U_{ijij} - U_{ijji} \delta_{\sigma\tau}
    """
    nbands = ufull.shape[0]
    udens = np.zeros((nbands, 2, nbands, 2), ufull.dtype)

    for b1 in range(0,nbands):
       for s1 in range(0,2):
          for b2 in range(0,nbands):
             for s2 in range(0,2):
                udens[b1,s1,b2,s2] += ufull[b1,s1,b2,s2,b1,s1,b2,s2]
                if s1 == s2:
                   udens[b1,s1,b2,s2] -= ufull[b1,s1,b2,s2,b2,s2,b1,s1]

    return udens


class Interaction:
   """Base class for storing the local interaction (the U matrix)."""
   def __init__(self,nbands):
      """Constructor of the base class."""
      self.norbitals = nbands
      self.u_matrix = np.zeros(shape=(nbands,2,nbands,2,nbands,2,nbands,2),
                               dtype=float,order='C')
      self.quantum_numbers = "Nt",

   def set_u_matrix(self):
      """Abstract Method to fill U-matrix with arbitrary interaction"""
      raise Exception('Abstract method, please override')

   def get_udens(self, Uw=0, Uw_Mat=0):
      """Abstract Method to get density-density part of interaction"""
      raise Exception('Abstract method, please override')

   def get_ubar(self):
       """Returns average U-value for use in double counting"""
       Uddav = 0.
       for im in range(self.norbitals):
           for imp in range(self.norbitals):
               Uddav += (self.u_matrix[im,0,imp,1,im,0,imp,1]+
                         self.u_matrix[im,0,imp,0,im,0,imp,0]-
                         self.u_matrix[im,0,imp,0,imp,0,im,0])
       Uddav /= self.norbitals * (2*self.norbitals-1)
       return Uddav

   def convert_u_matrix(self):
      """Method to convert the 8 dimensional python u-matrix to the four dimensional fortran u-matrix"""
      return self.u_matrix.reshape(self.norbitals*2,self.norbitals*2,self.norbitals*2,self.norbitals*2,order='C')

   def __eq__(self, other):
      """Abstract method to compare two interactions"""
      raise Exception('Abstract method, please override')

   def __ne__(self, other):
      return not self == other

   def dump(self,f=sys.stderr):
      """Dumping the U-matrix in the format required by readin."""

      print('# Non-zeros of interaction matrix U_ijkl', file=f)
      print(self.norbitals, 'BANDS', file=f)

      #extracting non-zero elements
      tmp = np.nonzero(self.u_matrix)
      pos = np.asarray(list(tmp)).T

      for i in pos:
         ind = tuple(i.tolist())
         ind_2 = ""
         for pos,val in enumerate(ind):
            #band index
            if pos%2==0:
               ind_2 += str(val+1)
            #spin index
            else:
               ind_2 += str(self._spin_char(val))
               ind_2 += " "

         print(ind_2, self.u_matrix[ind], file=f)

   def _char_spin(self, char):
      """Converting u,d spin to 0 or 1."""
      try:
         return {'u':0, 'd':1}[char]
      except KeyError:
         raise ValueError('Invalid spin, expecting u or d, got', char)

   def _spin_char(self, spin):
      """Converting 0,1 to u,d spin."""
      try:
         return ('u','d')[spin]
      except IndexError:
         raise ValueError('Invalid spin, expecting 0 or 1, got', spin)

class Density(Interaction):
   """Derived class for Density interaction."""
   def __init__(self, nbands, u, v, j):
      Interaction.__init__(self,nbands)
      self.name = "Density"
      self.u = u
      self.v = v
      self.j = j
      self.set_u_matrix()
      self.quantum_numbers = "Nt", "Szt", "Azt"

   def set_u_matrix(self):
      """Adding the density-density interaction to the U-matrix."""
      for b1 in range(self.norbitals):
         self.u_matrix[b1,(0,1),b1,(1,0),b1,(0,1),b1,(1,0)] = self.u
         for b3 in range(self.norbitals):
            if b3 != b1:
               self.u_matrix[b1,(0,0,1,1),b3,(0,1,0,1),b1,(0,0,1,1),b3,(0,1,0,1)] = self.v
               self.u_matrix[b1,(0,1),b3,(0,1),b3,(0,1),b1,(0,1)] = self.j

   def sym_u_matrix(self):
      """Adding crossing symmetric density-density U-matrix."""
      self.u_matrix=0.
      self.set_u_matrix()
      for b1 in range(self.norbitals):
         self.u_matrix[b1,(0,1),b1,(1,0),b1,(1,0),b1,(0,1)] += -self.u
         for b3 in range(self.norbitals):	
            if b3 != b1:
               self.u_matrix[b1,(0,0,1,1),b3,(0,1,0,1),b3,(0,1,0,1),b1,(0,0,1,1)] += -self.v
               self.u_matrix[b1,(0,1),b3,(0,1),b1,(0,1),b3,(0,1)] += -self.j
      
#      self.u_matrix*=0.5

   def __str__(self):
       return "<d-d: U=%g U'=%g J=%g>" % (self.u, self.v, self.j)

   def get_udens(self, Uw=0, Uw_Mat=0):
      if Uw==1:
        return screenedudensity_values(self.norbitals, Uw_Mat, self.u, self.v, self.j)
      else:
        return udensity_values(self.norbitals, self.u, self.v, self.j)

   def __eq__(self, other):
      """Comparing two U-Matrices by value."""
      if isinstance(other,self.__class__):
         return (self.u_matrix.shape[0] == other.u_matrix.shape[0] and self.u == other.u and
         self.v == other.v and self.j == other.j)
      else:
         return False

class Kanamori(Interaction):
    """Derived class for Kanamori interaction."""
    def __init__(self, nbands, u, v, j):
       Interaction.__init__(self, nbands)
       self.name = "Kanamori"
       self.u = u
       self.v = v
       self.j = j
       self.set_u_matrix()
       self.quantum_numbers = "Nt", "Szt", "Qzt"

    def set_u_matrix(self):
      """Adding the Kanamori interaction to the U-matrix."""

      self.u_matrix = Density(self.norbitals, self.u, self.v, self.j).u_matrix

      for b1 in range(self.norbitals):
         for b2 in range(self.norbitals):
            if b1 != b2:
               self.u_matrix[b1,(0,1),b1,(1,0),b2,(0,1),b2,(1,0)] = self.j # pair hopping
               self.u_matrix[b1,(0,1),b2,(1,0),b2,(0,1),b1,(1,0)] = self.j # spin flip

    def sym_u_matrix(self):
      """Adding crossing symmetric density-density U-matrix."""
      self.u_matrix=0.
      self.set_u_matrix()
      for b1 in range(self.norbitals):
         self.u_matrix[b1,(0,1),b1,(1,0),b1,(1,0),b1,(0,1)] += -self.u
         for b3 in range(self.norbitals):
            if b3 != b1:
               self.u_matrix[b1,(0,0,1,1),b3,(0,1,0,1),b3,(0,1,0,1),b1,(0,0,1,1)] += -self.v
               self.u_matrix[b1,(0,1),b3,(0,1),b1,(0,1),b3,(0,1)] += -self.j

               self.u_matrix[b1,(0,1),b1,(1,0),b3,(1,0),b3,(0,1)] += -self.j # pair hopping
               self.u_matrix[b1,(0,1),b3,(1,0),b1,(1,0),b3,(0,1)] += -self.j # spin flip
      
#      self.u_matrix*=0.5

       

    def __str__(self):
       return "<Kan: U=%g U'=%g J=%g>" % (self.u, self.v, self.j)

    def get_udens(self, Uw=0, Uw_Mat=0):
       return udensity_values(self.norbitals, self.u, self.v, self.j)

    def __eq__(self, other):
      """Comparing two U-Matrices by value."""
      if isinstance(other,self.__class__):
         return (self.u_matrix.shape[0] == other.u_matrix.shape[0] and self.u == other.u and
         self.v == other.v and self.j == other.j)
      else:
         return False

class Coulomb(Interaction):
   """Adding the Coulomb interaction to the U-matrix."""
   def __init__(self,nbands,f0,f2,f4,f6):
      Interaction.__init__(self,nbands)
      self.name = "Coulomb"
      self.f0 = f0
      self.f2 = f2
      self.f4 = f4
      self.f6 = f6
      self.set_u_matrix()
      self.quantum_numbers = "Nt", "Szt", "All"

   def Jplus(self,j,m):
      return np.sqrt(j*(j+1) - m*(m+1))
   def Jminus(self,j,m):
      return np.sqrt(j*(j+1) - m*(m-1))
     
   def CG(self,j1,m1,j2,m2,J,M):
      # this is an adaption from a code of Ragnar Stroberg
      if (m1+m2) != M or abs(m1)>j1 or abs(m2)>j2 or abs(M)>J:
          return 0
      if j1>j2+J or j1<abs(j2-J):
          return 0
      # Step 1:  Find all overlaps < (j1, k1)(j2,J-k1) | JJ >
      # They are related by applying J+ to the j1 and j2 components
      # and setting the sum equal to J+|JJ> = 0. Normalization is from completeness.
      # Start with < (j1,j1)(j2,J-j1) | JJ >, with weight 1
      n = 1
      Norm = [1]
      K = np.arange(j1,J-j2-1,-1)    
      for k in K[1:]:
          n *= -self.Jplus(j2,J-k-1) / self.Jplus(j1,k)
          Norm.append(n)       
      Norm /= np.sqrt(np.sum(np.array(Norm)**2))
      # Step 2: Apply J- to get from |JJ> to |JM>, and do the same
      # with j1 and j2. Do this for all the overlaps found in Step 1
      cg = 0
      for i,k1 in enumerate(K):
          k2 = J-k1
          if k1<m1 or k2<m2: continue
          # multiply by a factor F to account for all the ways
          # you can get there by applying the lowering operator
          F = factorial(k1-m1+k2-m2) / ( factorial(k1-m1) * factorial(k2-m2) )
          # Apply the lowering operators
          c1,c2,C = 1,1,1
          for k in np.arange(k1,m1,-1):
              c1 *= self.Jminus(j1,k)
          for k in np.arange(k2,m2,-1):
              c2 *= self.Jminus(j2,k)
          for k in np.arange(J,M,-1):
              C  *= self.Jminus(J,k)
          cg += c1*c2/C * Norm[i] * F       
      return cg
   def aFak(self,k, l, m1, m2, m3, m4):
       a = 0.

       for q in range(-k,k+1):
           a += (self.CG(l, m3, k, q, l, m1)*
                 self.CG(l, m2, k, q, l, m4))

       a = a * self.CG(l,0, k,0, l,0 )**2
       return a


   def set_u_matrix(self):
       """Initializing Coulomb interaction from the F0 ... F6 Slater parameters."""
       if np.sum(self.norbitals == np.array([1,3,5,7])) == 0:
           raise NotImplementedError('Coulomb only implemented for norbitals = 1, 3, 5, 7, i.e., s, p, d, f orbitals')
       l = (self.norbitals - 1)//2
       uuCub = np.zeros(shape=(self.norbitals,self.norbitals,self.norbitals,self.norbitals),dtype=complex)
       uuSph = np.zeros(shape=(self.norbitals,self.norbitals,self.norbitals,self.norbitals))
       mList = range(-l,l+1)
       for m1 in mList:
           for m2 in mList:
               for m3 in mList:
                   for m4 in mList:
                       # +l to fix the indexing   
                       uuSph[m1+l,m2+l,m3+l,m4+l]= (
                                self.f0 * self.aFak(0,l,m1,m2,m3,m4)
                              + self.f2 * self.aFak(2,l,m1,m2,m3,m4)
                              + self.f4 * self.aFak(4,l,m1,m2,m3,m4)
                              + self.f6 * self.aFak(6,l,m1,m2,m3,m4)
                                                )
       # transformation matrix from spherical harmonics to cubic harmonics

       rotMat = np.zeros(shape=(self.norbitals,self.norbitals),dtype=complex)
       sq2 = 1.0/np.sqrt(2)
       if self.norbitals == 7:
           # order of sphericals: m=-3,-2,-1,0,1,2,3
           # order of cubics: x(x^2-3y^2), z(x^2-y^2), xz^2, z^3, z^3, yz^2, xyz,  y(3x^2 - y^2)
           rotMat[0,0] = sq2
           rotMat[0,6] =-sq2
           rotMat[1,1] = sq2
           rotMat[1,5] = sq2
           rotMat[2,2] = sq2
           rotMat[2,4] =-sq2
           rotMat[3,3] = 1.0
           rotMat[4,2] = 1j*sq2
           rotMat[4,4] = 1j*sq2
           rotMat[5,1] = 1j*sq2
           rotMat[5,5] =-1j*sq2
           rotMat[6,0] = 1j*sq2
           rotMat[6,6] = 1j*sq2   
               
       elif self.norbitals == 5:
           # order of sphericals: m=-2, -1, 0, 1, 2
           # order of cubics: xy, yz, z^2, xz, x^2-y^2 (like in VASP)
           rotMat[0,0] = 1j*sq2;
           rotMat[0,4] =-1j*sq2
           rotMat[1,1] = 1j*sq2
           rotMat[1,3] = 1j*sq2
           rotMat[2,2] = 1.0
           rotMat[3,1] = sq2
           rotMat[3,3] =-sq2
           rotMat[4,0] = sq2
           rotMat[4,4] = sq2
       elif self.norbitals == 3:
           # order of sphericals: m=1,0,-1
           # order of cubics: y, z, x
           rotMat[0,0] = 1j*sq2
           rotMat[0,2] = 1j*sq2
           rotMat[1,1] = 1.0
           rotMat[2,0] = sq2
           rotMat[2,2] = -sq2
       elif self.norbitals == 1:
           rotMat[0,0] = 1                                              

       for b0 in range(self.norbitals):
           for b1 in range(self.norbitals):
               for b2 in range(self.norbitals):
                   for b3 in range(self.norbitals):
                      dum1 = np.tensordot(np.conj(rotMat[b0,:]),uuSph,axes=(0,0))
                      dum2 = np.tensordot(np.conj(rotMat[b1,:]),dum1,axes=(0,0))
                      dum3 = np.tensordot(dum2,rotMat[b2,:],axes=(0,0))
                      uuCub[b0,b1,b2,b3] = np.tensordot(dum3,rotMat[b3,:],axes=(0,0))
                      
       self.uuCub = uuCub   # save the spin-independent u-matrix              

       # make spin dependent u-matrix
       for b0 in range(self.norbitals):
          for b1 in range(self.norbitals):
              for b2 in range(self.norbitals):
                  for b3 in range(self.norbitals):
                      self.u_matrix[b0, (0,0,1,1), b1, (0,1,0,1), b2, (0,0,1,1), b3, (0,1,0,1)] = uuCub[b0,b1,b2,b3].real
  

   def __str__(self):
      return "<Coul: Fn=%g, %g, %g, %g>" % (self.f0, self.f2, self.f4, self.f6)

   def get_udens(self, Uw=0, Uw_Mat=0):
      return udensity_from_ufull(self.uuCub) 

   def __eq__(self, other):
      """Comparing two U-Matrices by value."""
      if isinstance(other,self.__class__):
         return (self.u_matrix.shape[0] == other.u_matrix.shape[0] and self.f0 == other.f0 and 
         self.f2 == other.f2 and self.f4 == other.f4 and self.f6 == other.f6)
      else:
         return False

class CustomFull(Interaction):
   """Read-In U Matrix from File class"""

   def __init__(self, u_matrix):
      Interaction.__init__(self, u_matrix.shape[0])
      self.u_matrix = u_matrix
      self.set_u_matrix()
      self.name = "CustomFull"
      self.quantum_numbers = "Nt", "All"

   def set_u_matrix(self):
      """Setting values of the U-matrix based on a file."""
      pass

   def get_udens(self, Uw=0, Uw_Mat=0):
      return udensity_from_spindep_ufull(self.u_matrix)

   def __eq__(self, other):
      if isinstance(other,self.__class__):
         return np.array_equal(self.u_matrix,other.u_matrix)
      else:
         return False

class CustomSU2Invariant(CustomFull):
   """Read-In U Matrix without spin index from File class"""
   def __init__(self, orbital_u_matrix):
      Interaction.__init__(self, orbital_u_matrix.shape[0])
      self.orbital_u_matrix = orbital_u_matrix
      self.set_u_matrix()
      self.name = "CustomSU2Invariant"
      self.quantum_numbers = "Nt", "Szt", "All"

   def set_u_matrix(self):
      """Setting values of the U-matrix based on a file."""
      # make spin-dependent one (honor SU(2) symmetry)
      self.u_matrix[:, 0, :, 0, :, 0, :, 0] = self.orbital_u_matrix
      self.u_matrix[:, 0, :, 1, :, 0, :, 1] = self.orbital_u_matrix
      self.u_matrix[:, 1, :, 0, :, 1, :, 0] = self.orbital_u_matrix
      self.u_matrix[:, 1, :, 1, :, 1, :, 1] = self.orbital_u_matrix

   def get_udens(self, Uw=0, Uw_Mat=0):
      return udensity_from_ufull(self.orbital_u_matrix)

def _close(a,b,rtol=1.e-05,atol=1e-16):
   """helper function to mimic np.allclose behavior for scalars."""
   return (abs(a-b) <= (atol + rtol*abs(b)) )

def similar(i1, i2, eps):
   """function to compare two interaction objects within eps
       for exact comparison refer to operator == != overloads in classes"""
   if isinstance(i1,i2.__class__):
      if isinstance(i1,Kanamori) or isinstance(i1,Density):
         return (i1.u_matrix.shape[0] == i2.u_matrix.shape[0] and
         _close(i1.u,i2.u,rtol=0.,atol=eps) and
         _close(i1.v,i2.v,rtol=0.,atol=eps) and
         _close(i1.j,i2.j,rtol=0.,atol=eps))
      if isinstance(i1,Coulomb):
         return (i1.u_matrix.shape[0] == i2.u_matrix.shape[0] and
         _close(i1.f0,i2.f0,rtol=0.,atol=eps) and
         _close(i1.f2,i2.f2,rtol=0.,atol=eps) and
         _close(i1.f4,i2.f4,rtol=0.,atol=eps) and
         _close(i1.f6,i2.f6,rtol=0.,atol=eps))
      if isinstance(i1,CustomFull) or isinstance(i1,CustomSU2Invariant):
         return np.allclose(i1.u_matrix,i2.u_matrix,rtol=0.,atol=eps)
   else:
      return False

