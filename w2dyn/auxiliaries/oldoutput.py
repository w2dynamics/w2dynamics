"""@package oldoutput
Provides class QmcOutput, which is used to convert old text output files to
a python object convenient for analysis in pylab.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import re
import warnings
import numpy as np
import configobj as cfo
from sys import stdout

class MatchFile(object):
    """ 
    Decorator of the file object for reading only lines matching some pattern.
    
    This object allows you to split "compound" files into their individual
    parts, if these parts are denoted by some special line format (e.g., a
    a hash sign as first character or having some content other than 
    whitespace).  The read(), readline() and readlines() methods will until they
    encounter a line not matching the pattern, in which case EOF is returned.
    You will then typically change the pattern or skip separation lines to read
    the next part of the file.
    
    Unlike the file object, this object allows you to mix iterator and read()
    API, but is rather slow compared to the native Python I/O API.  The class 
    follows the decorator pattern and exposes all methods of the underlying file
    class.  When using writing and seeking, though, be careful to leave the file
    pointer at the beginning of a line.
    """
    def __iter__(self):
        return self;
    
    def __init__(self, file, pattern=re.compile("")):
        """ Wraps a new MatchFile object over an existing file """
        if file.isatty():
            raise ValueError("File must be seek-able (no terminal)")
        self._file = file
        self.pattern = pattern
        
    def __next__(self):
        """ For use in the line iterators (for line in file:) """
        line = self.readline()
        if line == "":
            raise StopIteration()
        return line

    def next(self):
        return self.__next__()

    def peekline(self):
        """ Retrieves the next line without changing the file pointer """
        pos = self._file.tell();
        line = self._file.readline();
        self._file.seek(pos);
        return line;

    def readline(self, size=None):
        """ 
        Reads the next line if it matches the pattern, returns EOF otherwise.

        @param size: The maximum line size to be read.
        @raise ValueError:  
            If the line is longer than @c size, because this means that we
            would read an incomplete line in the next turn.  Before raising,
            the file pointer is reset to the beginning of the line to ensure
            a "clean" state.
        """
        pos = self._file.tell()
        if size is None:
            line = self._file.readline()
        else:
            line = self._file.readline(size)
        if line == "":
            return line;  # end of file
        if len(line) == size and not line.endswith("\n"):
            self._file.seek(pos) # revert to good state
            raise ValueError("Incomplete line read")
        # now do the match
        if self.pattern.match(line):
            return line
        else:
            self._file.seek(pos)
            return ""

    def readlines(self, sizehint=None):
        """ Returns all lines matching the pattern as a list """
        lines = []
        sizeread = 0
        while True:
            if sizehint and sizeread >= sizehint:
                break
            line = self.readline()
            if line == "":
                break
            lines.append(line)
        
        return lines;
        
    def read(self, size=None):
        """ Returns all matching lines as a string """
        try:
            from io import StringIO
        except ImportError:
            from StringIO import StringIO
        line = "s"
        output = StringIO()
        
        while line != "":
            if size is None:
                line = self.readline()
            else:
                line = self.readline(size);
                size = size - len(line);
                output.write(line); # also OK for EOF
            
        return output.getvalue()
    
    def __getattr__(self, name):
        return getattr(self._file, name)
#

class IterationData(object):
    """ Object holding the data of one CT-QMC iteration as extracted from
        the output file """
    pass

class Element(object):
    """ Data element """
    pass

class ArrayBuilder:
    """
    Compiles multi-dimensional arrays from a set of (index, value) pairs.
    """
    def __init__(self, ndims, elements, dimnames, nelements=None, elemtypes=None):
        self.ndims = ndims
        self.axisvals = [[] for i in range(ndims)]
        self.fixedsize = nelements is not None
        if self.fixedsize:
            assert elemtypes is not None, "elemtypes must be set if nelements is set"
            self._count = 0
            self._target = nelements
            self.values = dict( ((ename,np.zeros(shape=(nelements,),dtype=etype))
                                 for ename, etype in zip(elements,elemtypes)) )
        else:
            self.values = dict( ((elem,[]) for elem in elements) )
        self.dimnames = dimnames
        self._prev = None
        self._lockdim = ndims-1
    
    def append(self, indices, **values):
        """ Appends at the back of the array, honouring continuity """
        assert len(indices) == self.ndims, "Length not consistent"
        
        if self._prev == None:
            # initiates axis values
            self.axisvals = [ [index] for index in indices ]
        else:
            # Assertions guarantee right indexing
            assert indices != self._prev, "Duplicate index"
            dim,val = next(((i,v) for i,v in enumerate(indices) if v!=self._prev[i]))
            
            for i in range(dim+1, self.ndims):
                # check beginning of locked dim and end of old
                assert indices[i] == self.axisvals[i][0], "Illegal dimension opening"
                assert self._prev[i] == self.axisvals[i][-1], "Unclosed dimension"
            if dim > self._lockdim:
                # iterate in locked dimension
                assert self.axisvals[dim].index(val) - \
                       self.axisvals[dim].index(self._prev[dim]) == 1, "Non continuous"
                pass
            else:
                # extend unlocked dimension
                assert val not in self.axisvals[dim], "Skip-back indexing"
                self.axisvals[dim].append(val)
                self._lockdim = dim
        
        self._prev = indices
        if self.fixedsize:
            assert self._count < self._target, "Exceeded nelements"
            for k in self.values:
                self.values[k][self._count] = values[k]
            self._count += 1
        else:
            for k in self.values:
                self.values[k].append(values[k])
            
    def compile(self):
        """ Compiles the built-in arrays into numpy arrays """
        shape = [ len(axis) for axis in self.axisvals ]
        out = Element()

        for k in self.values:
            array = np.array(self.values[k]).reshape(shape)
            setattr(out, k, array)
        for name,vals in zip(self.dimnames, self.axisvals):
            setattr(out, name, vals)
        return out

_re_prefix = re.compile(":((?:[a-zA-Z_]|\d[a-zA-Z_])+)(\d*):")

def _chkfmt(name, pind, rind, pnum, rnum):
    """ Helper checking the right number of indices and values"""
    if pind > rind:
        warnings.warn("%s: Too many indices (%d > %d)" % (name, pind, rind))
    elif pind < rind:
        raise ValueError("%s: Too few indices  (%d < %d)" % (name, pind, rind))
    if pnum > rnum:
        warnings.warn("%s: Too many values (%d > %d)" % (name, pnum, rnum))
    elif pnum < rnum:
        raise ValueError("%s: Too few values  (%d < %d)" % (name, pnum, rnum))

def _tonum(str):
    try:
        return int(str)
    except ValueError:
        return float(str)

def _valueit(lineit):
    """ Wrapper routine over the line iterator extracting the indices from the
        prefix and the other values, returning them as a combined list """
    for i,line in enumerate(lineit):
        if i % 1000 == 999: 
           stdout.write(".")
           if i % 100000 == 99999:
               print()
           stdout.flush()
           
        toks = line.split()
        indices = [int(digit) for digit in _re_prefix.match(toks[0]).group(2)]
        numbers = [_tonum(tok) for tok in toks[1:]]
        yield indices + numbers
    if i >= 999:
        print()

class QmcOutput(IterationData):
    """ Structure for parsing and holding text output of the CT-QMC algorithm.
    
    The output file must have the following structure to be able to be parsed
    by the algorithm:
      1. It has a preamble containing set of key=value pairs (INI style), which 
         ends with the first line beginning with ":"
      2. It then continues with a list of lines with the following layout:
           <line>    ::= ":" <key> <indices> ":" <number>*
           <key>     ::= CHAR <key> | CHAR DIGIT* <key>
           <indices> ::= DIGIT*
      3. Consecutive lines with the same key are treated as single object,
         whereas multiple, non-consecutive instances of the same key are
         assigned to different iterations  
    
    Example usage:
    
       import OldOutput
       output = OldOutput.QmcOutput(open("output_file.out", "r"))
       print(output.data.GTau.mean[0,0,0,:])  # plot G(tau)
       import pickle
       pickle.dump(output, open("output_file.p","w"), -1)
    
    @param fname: Filename of the text input file.
    @returns: Object with three main attributes:
        - unknown -- keys not recognised
        - config -- the parsed text configuration in the file preamble
        - data -- the data of the final iteration, contains keys as attributes
        - it -- the data of all other iterations in the same fashion
    """
    def __init__(self, fname):
        #
        # Read configuration
        mf = MatchFile(open(fname,"r"), re.compile("^[^:]"))
        self.config = cfo.ConfigObj(infile=mf.readlines(), indent_type="\t")
        # 
        self.it = []
        self.unknown = []
        iter = IterationData()
        while True:
            # Look at the next line, and extract the prefix 
            line = mf.peekline()
            if line == "":
                break;
            match = _re_prefix.match(line)
            if not match:
                raise ValueError("Line has no prefix:\n"+line)

            id = match.group(1)
            digits = match.group(2)
            if hasattr(iter, id) and getattr(iter, id) is not None:
                # if we had that before, then start a new iteration
                self.it.append(iter)
                iter = IterationData()
                print("New iteration: ", len(self.it))
            
            # reset the pattern in order to read the right scope
            print("Found id:", id)
            mf.pattern = re.compile(":"+id+"[^a-zA-Z_]")
            value = None
            
            if id in ("NProcessors", "Gtotdensold", "Gtotdensnew", "mu", "QMCtotdens", "time"):
                # Handle scalar quantities
                for v in _valueit(mf):
                    assert value == None, "More than one value provided"
                    value = v
            
            elif id in ("mean_sign", "acc_add", "acc_rem", "acc_glob"):
                # Handle non-equivalent atom indexed quantities
                ab = ArrayBuilder(1, ['mean'], ['ineq'])
                for v in _valueit(mf):
                    ab.append(v[0:1], mean=v[1])
                value = ab.compile()
            
            elif id in ("Fiw", "Giw", "Siw"):
                # Handle non-equivalent atom indexed quantities
                ab = ArrayBuilder(4, ['mean'], ['ineq','band','spin','iwn'])
                for v in _valueit(mf):
                    ab.append((v[0],v[1],2*v[2]-3,v[3]), mean=complex(v[4],v[5]))
                value = ab.compile()
            
            elif id == "occ":
                ab = ArrayBuilder(5, ['mean','error'], 
                                  ['ineq','band1','spin1','band2','spin2'])
                for v in _valueit(mf):
                    ab.append((v[0], v[1], 2*v[2]-3, v[3], 2*v[4]-3), 
                              mean=v[5], error=v[6])
                value = ab.compile()
            
            elif id == "FTau":
                ab = ArrayBuilder(4, ['mean'], ['ineq','band','spin','tau'])
                for v in _valueit(mf):
                    ab.append((v[0], v[1], 2*v[2]-3, v[3]), mean=v[4])
                value = ab.compile()

            elif id in "GTau":
                ab = ArrayBuilder(4, ['mean','error'], ['ineq','band','spin','tau'])
                for v in _valueit(mf):
                    ab.append((v[0], v[1], 2*v[2]-3, v[3]), mean=v[4], error=v[5])
                value = ab.compile()
            
            elif id in "GTauLeg":
                ab = ArrayBuilder(4, ['mean'], ['ineq','band','spin','tau'])
                for v in _valueit(mf):
                    ab.append((v[0], v[1], 2*v[2]-3, v[3]), mean=v[4])
                value = ab.compile()
            
            elif id == "GLeg":
                ab = ArrayBuilder(4, ['mean','error'], ['ineq','band','spin','legorder'])
                for v in _valueit(mf):
                    ab.append((v[0], v[1], 2*v[2]-3, v[3]), mean=v[4], error=v[5])
                value = ab.compile()

            elif id == "Histo":
                ab = ArrayBuilder(4, ['share'], ['ineq','band','spin','traceorder'])
                for v in _valueit(mf):
                    ab.append((v[0], v[1], 2*v[2]-3, v[3]), share=v[4])
                value = ab.compile()
            
            elif id in ("GTau4Pnt", "GGTau", "Chi"):
                nbands = int(self.config["NBands"])
                ab = ArrayBuilder(8, ['mean','error'], ['ineq','band1','spin1','band2',
                                                        'spin2','tau12','tau34','tau14'])
                for v in _valueit(mf):
                    b1,b2 = (v[1]-1)//nbands+1, (v[1]-1)%nbands+1
                    s1,s2 = (v[2]-1)//2*2-1, (v[2]-1)%2*2-1
                    ab.append((v[0], b1, s1, b2, s2, v[3], v[4], v[5]), mean=v[6], error=v[7])
                value = ab.compile()
            
            elif id == "GLeg4Pnt":
                nbands = int(self.config["NBands"])
                ab = ArrayBuilder(8, ['mean','error'], ['ineq','band1','spin1','band2',
                                                        'spin2','l','lprime','iw'])
                for v in _valueit(mf):
                    b1,b2 = (v[1]-1)//nbands+1, (v[1]-1)%nbands+1
                    s1,s2 = (v[2]-1)//2*2-1, (v[2]-1)%2*2-1
                    ab.append((v[0], b1, s1, b2, s2, v[3], v[4], v[5]), mean=complex(v[6],v[7]),
                              error=v[8])
                value = ab.compile()
            
            # Handle unrecognized entries
            else:
                mf.readlines()
                if id not in self.unknown:
                    self.unknown.append(id)
                    warnings.warn("Unrecognized entry: " + id)
            
            setattr(iter, id, value)
            #iter.id.__desc__ = "No description available"
        # 
        mf.close()
        self.it.append(iter)
        self.data = iter
        
        try:
            self.beta = float(self.config["beta"])
            self.nbands = int(self.config["NBands"])
            self.niw = int(self.config["NGiw"])
            self.ntau = int(self.config["NGtau"])
            self.legorder = int(self.config["NLegOrder"]) 
        except KeyError:
            warnings.warn("Did not find common key in config")
            
