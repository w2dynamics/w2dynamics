"""This module contains conversion functions for band-spin patterns into 
compound indices and vice-versa.

Calculations are done in the class GFComponent.

The subsequent wrapper functions mirror the Fortran subroutine defined
CompoundIndex.F90.

Authors: Josef Kaufmann, Patrik Gunacker""" 
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np


# class for Green's function component
class GFComponent(object):
    """ Class for indexing green's functions.
    An instance of GFComponent holds the following fields:
    *** index: compound index of all indices (one-based single number)
    *** bands: band indices (zero-based list)
    *** spins: spin indices (zero-based list)
    *** bandspin: band-spin compound indices
    *** n_bands: number of impurity orbitals
    *** n_ops: number of operators in Greens function 
             2 -> 1-particle Green's function
             4 -> 2-particle Green's function"""


    def __init__(self, index=None, 
                 bands=None, spins=None, bandspin=None,
                 n_bands=0, n_ops=0, n_spins=2):

        if n_bands == 0 or n_ops == 0:
            raise ValueError('n_bands and n_ops have to be set'
                           + ' to non-zero positive integers')

        self.n_bands = n_bands
        self.n_ops = n_ops
        dims_bs = n_ops * (n_bands*n_spins,)
        dims_1 = (n_bands, n_spins)

        if index is not None and bands is None: # initialize from compound index
            self.index = index
            self.bandspin = list(np.unravel_index(self.index-1, dims_bs))
            self.bands,self.spins = np.unravel_index(self.bandspin, dims_1)

        elif bands is not None and index is None: # initialize from bands (and spins)
            self.bands = bands
            if spins is None: # use only band indices (e.g. d/m channel)
                self.spins = n_ops * (0,)
            else:
               self.spins = spins

            self.bandspin = np.ravel_multi_index(
                (self.bands,self.spins), (n_bands,n_spins))
            self.index = np.ravel_multi_index(self.bandspin, dims_bs) + 1

        elif bandspin is not None and index is None:
            self.index = np.ravel_multi_index(bandspin, dims_bs) + 1

        else:
            raise ValueError('index and bands both supplied')


    def bsbs(self):
        bsbs = np.vstack((self.bands,self.spins)).transpose().reshape(-1)
        return tuple(bsbs)


def index2component_general(Nbands, N, ind):
    """ converting an index into a band-spin pattern
    :param N: number of operators
    :param ind: general flavor index
    :return bandspin: band-spin array of length N
    :return bands: band array of length N
    :return spins: spin array of length N"""
    
    comp = GFComponent(index=ind, n_bands=Nbands, n_ops=N)
   
    return comp.bandspin, comp.bands, comp.spins

def index2component_band(Nbands, N, ind):
    """ converting an index into a band pattern
    :param N: number of operators
    :param ind: general flavor index
    :return bands: band array of length N"""
    
    comp = GFComponent(index=ind, n_bands=Nbands, n_ops=N, n_spins=1)
   
    return comp.bands

def component2index_general(Nbands, N, b, s):
    """ converting a band-spin pattern into an index
    :param N: number of operators
    :param b: band array of length N
    :param s: spin array of length N
    :return index: general flavor index"""
      
    comp = GFComponent(n_bands=Nbands, n_ops=N, bands=b, spins=s)
     
    return comp.index

def componentBS2index_general(Nbands, N, bs):
    """ converting a bandspin pattern into an index
    :param N: number of operators
    :param bs: bandspin array of length N
    :return index: general flavor index"""
    
    comp = GFComponent(n_bands=Nbands, n_ops=N, bandspin=bs)
     
    return comp.index

def component2index_band(Nbands, N, b):
    """ converting a band pattern into an index
    :param N: number of operators
    :param b: band array of length N
    :return index: general flavor index"""

    comp = GFComponent(n_bands=Nbands, n_ops=N, n_spins=1, bands=b)
     
    return comp.index
