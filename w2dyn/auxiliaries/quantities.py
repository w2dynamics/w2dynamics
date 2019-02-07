"""Tools for selecting data from HDF output files"""
import sys
import re
import os
import itertools
import numpy as np
import h5py as hdf5
from warnings import warn

debug = None

try:
    from compatibility import iter_product
    from postprocessing import derived_quantities
except ImportError:
    warn("import problems, reduced functionality")
    from itertools import product as iter_product
    derived_quantities = {}

def ineq_quantity(iter, qname, field=("value", "error"), ineq=None, 
                  file_version=None):
    if file_version is None:
        file_version = tuple(iter.parent.attrs["outfile-version"])
    try:
        ineq = np.atleast_1d(np.asarray(ineq, dtype=np.int))
    except TypeError:
        if ineq is None: ineq = slice(None)
    if qname in iter:
        qnode = iter[qname]
        nineqs = qnode["value"].shape[0]
        for iineq in np.arange(nineqs)[ineq]:
            yield tuple(qnode[key][iineq] if key in qnode else None
                        for key in field)
    elif file_version >= (2,0):
        key_pool = np.array(sorted(k for k in iter.keys() 
                                   if k.startswith("ineq-")), dtype=np.string_)
        for ineq_key in key_pool[ineq]:
            qnode = iter[ineq_key][qname]
            yield tuple(qnode[key].value if key in qnode else None
                        for key in field)
    else:
        raise ValueError("quantity not found in iteration")

# def get_config(hdf_file, cnames, ineq=None, nneq=None, file_version=None):
#     parse_value = Selector.parse_value
#     if file_version is None:
#         file_version = tuple(hdf_file.attrs["outfile-version"])
#     if nneq is None:
#         nneq = int(hdf_file["start/nneq"].value)
#     if file_version < (2,0):
#         cfg_dict = dict(PairListSelection._pairs(hdf_file["config"].value))
#         for cname in cnames:
#             qualifier, name = cname.rsplit(".", 1)
#             try:
#                 value = cfg_node.get(name)
#             except KeyError:
#                 yield None
#                 continue
#             value = parse_value(value)
#             if qualifier == "atoms.*":
#                 yield [value]*nneq
#             else:
#                 yield value
#     else:
#         cfg_dict = dict(hdf_file[".config"].attrs)
#         for cname in cnames:
#             qualifier, name = cname.rsplit(".", 1)

class Selector:
    """Essentially a tuple of single values or (extended) slices
    
    Schematics of selection:                               _____________ 
    ------------------------             _________________|             |
             ______________             |          _______|  Container  |     
            |              |            |         |   list|_____________|
            |   Selector   |_________   |         |         _____|___________       
            |______________|         |  |         |   _____|_______     _____|_____ 
             ______|______        ___V__V______   |  |   Wrapper   |   |           |
         ____|____    ____|____  |             |  ^  |  _________  |   | list/dict |
        |         |  |         | |  Selection  |  |  | |  .inner |---->|___________|
        | scalar  |  |  slice  | |_____________|  |  | |_________| |                
        |_________|  |_________|      |.items()   |  |_____________|                
                                      |___________|  
    """
    def __init__(self, *slices):
        self._slices = list(slices)
        
    def __str__(self):
        return ",".join(map(self._slice_to_string, self._slices)) 
    
    def __iter__(self):
        return self._slices.__iter__()
    
    def append(self, val):
        self._slices.append(val)

    def _slice_to_string(self, sl):
        def handle_none(val):
            if val is None: return ""
            return str(val)
        if isinstance(sl, slice):
            if sl.step is None:
                tup = sl.start, sl.stop
            else:
                tup = sl.start, sl.stop, sl.step
            return ":".join(map(handle_none, tup))
        else:
            return handle_none(sl)
    
    def transform(self, tf_value):
        """Applies a value transformation to a selector, returning the result"""
        tf = Selector()
        for sel in self._slices:
            if isinstance(sel, slice):
                tf.append(slice(tf_value(sel.start,0), tf_value(sel.stop,1),
                                tf_value(sel.step,2)))
            else:
                tf.append(tf_value(sel,0))
        return tf
    
    def empty(self):
        """Returns whether the selector contains slices"""
        return not self._slices
    
    @classmethod
    def parse_value(cls, val):
        """Helper routine for parsing scalar values for from_string"""
        val = val.strip()
        if val == "": return None
        try: return int(val)
        except ValueError:
            try: return float(val)
            except ValueError: return val
    
    @classmethod
    def from_string(cls, sel_string):
        """Parses a selector from a string, returning name and selector"""
        try:
            name, sel_string = sel_string.split("=", 1)
        except ValueError:
            name = None
        slices = []
        for sl in sel_string.split(","):
            parts = sl.split(":")
            if len(parts) == 1:
                val = cls.parse_value(parts[0])
                if val is None: raise ValueError("empty value")
                slices.append(val)
            else:
                slices.append(slice(*map(cls.parse_value, parts)))
        return name, cls(*slices)

class SelectorPool:
    """Pool of selectors that allows handling of positional arguments"""
    def __init__(self, *positional, **named):
        self.positional = list(positional)
        self.named = dict(named)
        
    def expect(self, key):
        """Fetches a selector by name, or falls back to next positional""" 
        try:
            return self.named[key]
        except KeyError:
            # HACK: this is because otherwise the break with old interface
            #       would be too severe
            if key == "field": return Selector()
            try:
                sel = self.positional.pop(0)
            except IndexError:
                sel = Selector()
            self.named[key] = sel
            return sel
        
    def add(self, key, value):
        """Adds a selector, with an optional key"""
        if key is None:
            self.positional.append(value)
        else:
            warn("duplicate selector: key already in selector stack")
            self.named[key] = value
    
    def __str__(self):
        return " ".join(map(str, self.positional) +
                        ["%s=%s" % (k,v) for k,v in self.named.iteritems()])

    def append_from_string(self, str, override=False):
        """Appends a selector by string"""
        name, sel = Selector.from_string(str)
        if name is None:
            self.positional.append(sel)
        else:
            if not override and name is self.named:
                raise RuntimeError("duplicate named selector %s" % name)
            self.named[name] = sel

class Selection:
    """Holds the result of a selector being applied to some container."""
    def __str__(self):
        return self.__class__.__name__
    def keys(self):
        """Returns the keys to the selections, i.e., the indices"""
        raise NotImplementedError()
    def values(self):
        """Returns the values to the selections, i.e.,. the data"""
        raise NotImplementedError()
    def items(self):
        """Returns enumeration of pairs, equivalent to zip(keys(), values())"""
        raise NotImplementedError()

class IndexedSelection(Selection):
    """Selection of a list of pairs by its index """
    def _select(self, list):
        if self.selector.empty(): return list
        selection = []
        for sel in self.selector:
            if isinstance(sel, slice):
                selection.extend(list[sel])
            else:
                selection.append(list[sel])
        return tuple(selection)
    
    def __init__(self, selector, list):
        self.selector = selector
        self._pairs = self._select(tuple(enumerate(list)))
    
    def __str__(self):
        return "%s(%s): %s" % (self.__class__.__name__[:-9], self.selector, 
                               self._pairs) 
        
    def keys(self):
        return tuple(pair[0] for pair in self._pairs)
    def values(self):
        return tuple(pair[1] for pair in self._pairs)
    def items(self):
        return self._pairs

class OneBasedSelection(IndexedSelection):
    """Inclusive one-based selection mimicking the behaviour of Fortran"""
    @classmethod
    def tf_value(cls, val, pos):
        """Transforms one-based value to zero-based one"""
        if val is None: return None
        val = int(val)
        if val == 0: raise ValueError("invalid zero in one-based indexing")
        if pos == 0:
            return val - int(val > 0)
        elif pos == 1:
            if val == -1: return None
            return val + int(val < 0)
        else:
            return val
        
    def __init__(self, selector, list):
        self.selector = selector.transform(OneBasedSelection.tf_value)
        self._pairs = self._select(tuple(enumerate(list, 1)))
    
class AxisSelection(IndexedSelection):
    """Selection on an ordered numeric axis"""
    @classmethod
    def transform(cls, selector, axis):
        def tf_val(val):
            if val is None: return None
            return np.searchsorted(axis, val)
        tf = Selector()
        for sel in selector:
            if isinstance(sel, slice):
                stride = sel.step
                if stride is not None: 
                    stride = int(np.round(stride/(axis[1]-axis[0])))
                tf.append(slice(tf_val(sel.start), tf_val(sel.stop), stride))
            else:
                tf.append(tf_val(sel))
        return tf

    def __init__(self, selector, axis, lst):
        self.selector = AxisSelection.transform(selector, axis)
        self._pairs = self._select(zip(axis, lst))

class IterationSelection(OneBasedSelection):
    """Selection of an iteration in a HDF5 file"""
    ITERRE = re.compile(r"(?:dmft|stat|worm)-\d+$")

    @classmethod
    def get_set(cls, f):
        try:
            iters = [("start", f["start"])]
        except KeyError:
            raise RuntimeError("output file has no start iteration")
        iters.extend(p for p in sorted(f.items()) if cls.ITERRE.match(p[0]))
        try:
            iters.append(("finish", f["finish"]))
        except KeyError:
            try:
                iters.append(("aborted", f["aborted"]))
            except KeyError:
                iters.append(("------", {"desc":""}))
        return iters

    @classmethod
    def tf_value(cls, val, pos):
        """Transforms iteration value to index in iteration array"""
        try:
            # parse string
            val = val.lower()
            if val == "start": return 0
            if val == "finish": return -1
            if val == "aborted": return -1
        except AttributeError:
            # squeeze interaction
            val = OneBasedSelection.tf_value(val, pos)
            if val is None: return (1, -1, None)[pos]
            return val + (-1,1)[val >= 0]

    @classmethod
    def tf_selector(cls, selector):
        tf = Selector()
        for sel in selector:
            if sel == slice(None):
                tf.append(sel)
            elif isinstance(sel, slice):
                tf.append(slice(cls.tf_value(sel.start,0), cls.tf_value(sel.stop,1),
                                cls.tf_value(sel.step,2)))
            else:
                tf.append(cls.tf_value(sel,0))
        return tf

    def __init__(self, selector, parent):
        assert parent.name == "/", "must be in root"
        self.selector = IterationSelection.tf_selector(selector)
        self._iters = IterationSelection.get_set(parent)
        self._pairs = self._select(self._iters)

class IneqSelection(OneBasedSelection):
    """Selection of an inequivalent atom group in the new file format"""
    INEQRE = re.compile(r"ineq-\d+")
    
    def __init__(self, selector, parent):
        self.selector = selector.transform(OneBasedSelection.tf_value)
        self._ineqs = tuple(p[1] for p in sorted(parent.items()) 
                            if IneqSelection.INEQRE.match(p[0]))
        self._pairs = self._select(tuple(enumerate(self._ineqs, 1)))

class DatasetSelectionFactory:
    """Prepares selections on arbitrary datasets of an agreeing shape"""
    @classmethod
    def _tf_sinkpos(cls, selector, dimsize):
        """Returns the corresponding sink slots for a selection""" 
        selection = Selector()
        pos = 0
        arr = range(dimsize)
        for sel in selector:
            if isinstance(sel, slice):
                sllen = len(arr[sel])
                selection.append(slice(pos,pos+sllen))
                pos += sllen
            else:
                selection.append(pos)
                pos += 1
        return selection, pos
    
    def __init__(self, dataset_shape, axes, *selectors):
        self.dataset_shape = dataset_shape
        ndim = len(dataset_shape)
        if axes is None: axes = [None]*ndim
        if len(axes) != ndim:
            raise ValueError("no. of axes (%d) vs. dataset shape (%s)" % 
                             (len(axes), dataset_shape))
        if len(selectors) > ndim:
            raise ValueError("too many selectors")
        if len(selectors) < ndim:
            if debug: debug("padding selectors with blancos")
            selectors = (list(selectors) +
                         [Selector()]*(ndim - len(selectors)))
        self.source_selectors = []
        self.sink_selectors = []
        self.sink_shape = []
        self.sink_axes = []
        for dimsize, selector, axis in zip(dataset_shape, selectors, axes):
            if selector.empty():
                selector = Selector(slice(None))
            elif axis is not None:
                selector = AxisSelection.transform(selector, axis)
            else:
                selector = selector.transform(OneBasedSelection.tf_value)
            if axis is None:
                axis = np.arange(1, dimsize+1)
            
            sink_sel, sink_dim = DatasetSelectionFactory._tf_sinkpos(selector, 
                                                                     dimsize)
            self.source_selectors.append(selector)
            self.sink_selectors.append(sink_sel)
            self.sink_shape.append(sink_dim)
            self.sink_axes.append(IndexedSelection(selector, axis).values())
        
        if debug: debug("sink_shape = %s", self.sink_shape)
        if debug: debug("sink_selectors = %s", map(str, self.sink_selectors))
        if debug: debug("sink_axes = %s", ["%s..." % str(a[:10]) 
                                           for a in self.sink_axes])
        
    def select(self, dataset):
        if dataset.shape != self.dataset_shape:
            raise ValueError("shape mismatch (wrong factory)")
        if not dataset.shape:
            sink = dataset[...]
        sink = np.empty(shape=self.sink_shape, dtype=dataset.dtype)
        if debug: sink[...] = -47119999
        for ssel, dsel in zip(iter_product(*self.sink_selectors), 
                              iter_product(*self.source_selectors)):
            sink[ssel] = dataset[dsel]
        return DatasetSelection(self.sink_axes, sink)

class DatasetSelection(Selection):
    """Generated by DatasetSelectionFactory for a concrete dataset"""
    def __init__(self, axes, data):
        self.axes = axes
        self.data = data
    def keys(self):
        return iter_product(*self.axes)
    def values(self):
        return self.data

class FileSelection(Selection):
    """Selection of HDF5 files"""
    @classmethod
    def get_latest(cls, directory='.', ext='.hdf5'):
        files = [f for f in os.listdir(directory) if f.endswith(ext)]
        if not files:
            raise ValueError("no HDF5 file in current directory")   
        return max(files, key=lambda f: os.stat(f).st_mtime)

    def __init__(self, selector, directory='.'):
        if selector.empty():
            raise ValueError("expected at least one filename or 'latest'")
        self.selector = []
        for f in selector:
            if isinstance(f, slice):
                raise ValueError("hgrep does not support ':' in hdf_file names")
            if f == 'latest':
                f = FileSelection.get_latest(directory)
                if debug: debug("latest hdf_file: %s", f)
            self.selector.append(f)
        
    def values(self):
        for f in self.selector:
            if debug: debug("next file: %s", f)
            f = hdf5.File(str(f), "r")
            yield f
            f.close()

class DictSelection(IndexedSelection):
    def __init__(self, selector, dct):
        self._pairs = []
        self.selector = selector
        if selector.empty():
            self._pairs = sorted(dct.items())
        for sel in selector:
            if isinstance(sel, slice):
                raise ValueError("ranges are not supported in named selection")
            sel = str(sel)
            if sel.find("*") >= 0:
                import fnmatch
                for k, v in dct.iteritems():
                    if fnmatch.fnmatch(k, sel): self._pairs.append((k,v))
            else:
                self._pairs.append((sel, dct.get(sel)))

class PairListSelection(DictSelection):
    @classmethod
    def _pairs(cls, lst):
        for line in lst:
            try:
                line = line[:line.index("#")]
            except ValueError:
                pass
            line = line.strip()
            if line:
                k, v = line.split("=", 1)
                yield k.rstrip(), Selector.parse_value(v)

    def __init__(self, selector, lst):
        DictSelection.__init__(self, selector, dict(PairListSelection._pairs(lst)))

class QuantityContainer:
    def _get_metadata(self):
        # get quantity metadata
        try:
            if self.file_version >= (2,0):
                if debug: debug("fetching quantity metadata (new way)")
                attrs = self.hdf_file[".quantities"][self.qname].attrs
            else:
                if debug: debug("fetching quantity metadata (old way)")
                attrs = next(group for key, group in 
                             IterationSelection.get_set(self.hdf_file)
                             if self.qname in group)[self.qname].attrs
        except (KeyError, StopIteration):
            raise KeyError("quantity `%s' not found in HDF5 file `%s'" % 
                           (self.qname, self.hdf_file.filename))
        self.meta = dict(attrs)
        if debug: debug("quantity metadata: %s", self.meta)
        # get axes metadata
        try:
            self.axes_names = list(attrs["axes"])
        except KeyError:
            self.axes_names = ()
        
    def _get_byineq(self):
        if (self.file_version >= (2,0) and self.axes_names and 
            self.axes_names[0] == "ineq"):
            if debug: debug("using by-ineq selection of new file format ")
            self.by_ineq = True
            self._prefix = (" iter", "ineq")
            self.axes_names.pop(0)
        else:
            self.by_ineq = False
            self._prefix = (" iter",)
        
    def _get_axes(self):
        self.axes_values = []
        axes_group = self.hdf_file[("axes", ".axes")[self.file_version[0]-1]]
        for axis_name in self.axes_names:
            try:
                self.axes_values.append(axes_group[axis_name].value)
                if debug: debug("found axis %s: %s...", axis_name, self.axes_values[-1][:10])
            except KeyError:
                self.axes_values.append(None)
                if debug: debug("did not find %s axis values", axis_name)
                
    def __init__(self, quantity_name, hdf_file, file_version, get_axes=True):
        self.qname = quantity_name
        self.hdf_file = hdf_file
        self.file_version = file_version
        self._get_metadata()
        if get_axes:
            self._get_byineq() 
            self._get_axes()

    def _return_qnode(self, parent):
        return parent[self.qname]

    def select(self, pool):
        iter_sel =  IterationSelection(pool.expect("iter"), self.hdf_file)
        for iter_name, iter_node in iter_sel.items():
            if debug: debug("entering iteration node: %s", iter_name)
            if self.by_ineq:
                ineq_sel = IneqSelection(pool.expect("ineq"), iter_node)
                for ineq_name, ineq_node in ineq_sel.items():
                    if debug: debug("entering ineq node: %s", ineq_name)
                    try:
                        yield (self._prefix, (iter_name, ineq_name), 
                               self._return_qnode(ineq_node))
                    except KeyError:
                        if debug: debug("did not find quantity %s", self.qname)
                        continue
            else:
                if debug: debug("searching directly in iteration node")
                try:
                    yield self._prefix, (iter_name,), self._return_qnode(iter_node)
                except KeyError:
                    if debug: debug("did not find quantity %s", self.qname)
                    continue 

class DerivedQttyContainer(QuantityContainer):
    def __init__(self, meta, hdf_file, file_version):
        if debug: debug("metadata: %s", meta)
        self.meta = meta
        self.hdf_file = hdf_file
        self.axes_names = meta["axes"]
        self.file_version = file_version
        self._get_byineq()
        self._get_axes()
        self.func = meta["func"]
        self.fields = meta["fields"]
        
        # base containers
        self.qname = "derived"
        self.base = meta["base"]
        
    def _return_qnode(self, parent):
        return tuple(parent[baseq] for baseq in self.base)
    
    def select(self, pool):
        for prefix, path, bqnodes in QuantityContainer.select(self, pool):
            if debug: debug("found derived bases %s", [n.name for n in bqnodes])
            valargs = []
            errargs = []
            for node in bqnodes:
                valargs.append(node["value"].value)
                try:
                    errargs.append(node["error"].value)
                except KeyError:
                    errargs.append(None)
            if debug: debug("arguments: %s", str(valargs + errargs)[:200])
            result = self.func(*(valargs + errargs))
            if not isinstance(result, tuple):
                result = (result,)
            if debug: debug("result: %s", str(result)[:200])
            yield prefix, path, dict(zip(self.fields, result))

class MetaQttyContainer:
    DESC = {
        "config": ("configfile",  ".config",      "Parameters of the run"),
        "env":    ("environment", ".environment", "shell environment of the run"), 
        "error":  ("error",       ".error",       "ERROR MESSAGE (python) of the run"),
        }

    def __init__(self, quantity_name, hdf_file, file_version):
        self.qname = quantity_name
        self.hdf_file = hdf_file
        self.file_version = file_version

        meta = MetaQttyContainer.DESC[self.qname]
        self.nodename = meta[self.file_version[0]-1]
        self.desc = meta[2]

    def select(self, pool):
        if self.file_version >= (2,0):
            return DictSelection(pool.expect("item"), 
                                 self.hdf_file[self.nodename].attrs)
        else:
            return PairListSelection(pool.expect("item"), 
                                     self.hdf_file[self.nodename].value)

class FieldContainer:
    def __init__(self, node):
        self.node = node
        self.shape = self.node.shape
        self.dtype = self.node.dtype
    def __getitem__(self, selector):
        return self.node[selector]

class RealContainer(FieldContainer):
    def __init__(self, node):
        FieldContainer.__init__(self, node)
        self.dtype = np.float
    def __getitem__(self, selector):
        return FieldContainer.__getitem__(self,selector).real

class ImagContainer(FieldContainer):
    def __init__(self, node):
        FieldContainer.__init__(self, node)
        self.dtype = np.float
    def __getitem__(self, selector):
        return FieldContainer.__getitem__(self,selector).imag

class FieldSelection(DictSelection):
    @classmethod
    def field_order(cls, field_tuple):
        try:
            return {"value":    0,
                    "value-re": 1,
                    "value-im": 2,
                    "spin-up":  30,
                    "spin-dn":  31,
                    "error":    100,
                   }[field_tuple[0]]
        except KeyError:
            if debug: debug("unknown field name")
            return 4711

    def __init__(self, selector, node):
        dct = {}
        self.shape = None
        for field_name, field_node in node.iteritems():
            self.shape = field_node.shape
            if issubclass(field_node.dtype.type, np.complexfloating):
                dct[field_name + "-re"] = RealContainer(field_node)
                dct[field_name + "-im"] = ImagContainer(field_node)
            else:
                dct[field_name] = FieldContainer(field_node)
        DictSelection.__init__(self, selector, dct)
        self._pairs = sorted(self._pairs, key=FieldSelection.field_order)
