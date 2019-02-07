import time
import os
from warnings import warn
import h5py as hdf5
import numpy as np

import w2dyn.auxiliaries as _aux
import w2dyn.auxiliaries.transform as _tf
import w2dyn.dmft.meta as meta
import w2dyn.dmft.orbspin as orbspin

def _flat_items(d, prefix=""):
    items = []
    for k, v in d.iteritems():
        fullk = prefix + k
        if isinstance(v, dict):
            items.extend(_flat_items(v, fullk + "."))
        else:
            items.append((fullk, v))
    return items

class HdfOutput:
    def load_old(self, oldfile, olditer=-1):
        hf = hdf5.File(oldfile, "r")
        file_version = tuple(hf.attrs["outfile-version"])

        if olditer == -1:
            oldit = hf["dmft-last"]
        else:
            oldit = hf["dmft-%03d" % olditer]

        ineqs = [node for name, node in sorted(oldit.iteritems())
                 if name.startswith("ineq-")]
        if not ineqs:
            raise ValueError("No impurity results found in iteration")

        mu = oldit["mu/value"].value
        smom_dd = None
        dc_latt = None

        if file_version >= (2, 1):
            siw_dd = [ineq["siw-full/value"].value.transpose(4,0,1,2,3)
                      for ineq in ineqs]
            smom_dd = [ineq["smom-full/value"].value.transpose(4,0,1,2,3)
                       for ineq in ineqs]
            try:
                dc_latt = oldit["dc-latt/value"].value
            except KeyError:
                pass
        elif file_version >= (2, 0):
            # FIXME: double-counting for old runs
            siw_dd = [orbspin.promote_diagonal(ineq["siw/value"].value.transpose(2,0,1))
                      for ineq in ineqs]
        else:
            # FIXME: double-counting for very old runs
            siw_dd = [orbspin.promote_diagonal(siw_ineq.transpose(2,0,1))
                      for siw_ineq in oldit["siw/value"].value]

        hf.close()
        return mu, siw_dd, smom_dd, dc_latt

    def __init__(self, config, git_revision=None, start_time=None, 
                 ineq_list=None, mpi_comm=None):
        # generate name for HDF5 file
        if start_time is None: start_time = time.localtime()
        runstring = time.strftime("%Y-%m-%d-%a-%H-%M-%S", start_time)
        self.filename = "%s-%s.hdf5" % (config["General"]["FileNamePrefix"],
                                        runstring)

        if mpi_comm is not None:
            self.is_writer = mpi_comm.Get_rank() == 0
        else:
            self.is_writer = True

        if self.is_writer:
            self.file = hdf5.File(self.filename, "w-")
            self.iter = self._ready_file(self.file, config, git_revision,
                                         start_time)
        else:
            self.file = None
            self.iter = None

        # store handles
        self.ineq_list = ineq_list
    
    def _ready_file(self, hf, config, git_revision, start_time):
        # generate HDF5 file and write important metadata
        hf.attrs["outfile-version"] = _aux.OUTPUT_VERSION
        hf.attrs["code-version"] = _aux.CODE_VERSION
        if git_revision is not None:
            hf.attrs["git-revision"] = git_revision
        hf.attrs["run-date"] = time.strftime("%c", start_time)

        # write configureation of current run
        hfgrp = hf.create_group(".config")
        for key, val in _flat_items(config):
            if val is None: continue
            hfgrp.attrs[key.lower()] = val

        # write environment of current run
        hfgrp = hf.create_group(".environment")
        for k,v in os.environ.iteritems():
            if v is None: continue
            hfgrp.attrs[k] = v

        # write metadata for possible axes objects
        ax = hf.create_group(".axes")
        beta = config["General"]["beta"]
        qcfg = config["QMC"]
        ax.create_dataset("iw", data=_tf.matfreq(beta, "fermi", 2*qcfg["Niw"]))
        ax.create_dataset("tau", data=np.linspace(0, beta, qcfg["Ntau"]))
        ax.create_dataset("tauf", data=np.linspace(0, beta, qcfg["Nftau"]))
        ax.create_dataset("taubin", data=_tf.tau_bins(beta, qcfg["Ntau"], 'centre'))
        ax.create_dataset("tausus", data=np.linspace(0, beta, qcfg["Ntau"]+1))
        ax.create_dataset("tau-g4", data=_tf.tau_bins(beta, qcfg["N4tau"], 'centre'))
        ax.create_dataset("iwb-g4", data=_tf.matfreq(beta, 'bose', 2*qcfg["N4iwb"]+1))
        ax.create_dataset("iwf-g4", data=_tf.matfreq(beta, 'fermi', 2*qcfg["N4iwf"]))
        ax.create_dataset("iwf-g2", data=_tf.matfreq(beta, 'fermi',
                                                     2*(qcfg["N4iwf"]+qcfg["N4iwb"])))
        ax.create_dataset("iwb-p2", data=_tf.matfreq(beta, 'bose',2*qcfg["N2iwb"]+1))
        ax.create_dataset("iwf-p3", data=_tf.matfreq(beta, 'fermi',2*qcfg["N3iwf"]))
        ax.create_dataset("iwb-p3", data=_tf.matfreq(beta, 'bose', 2*qcfg["N3iwf"]+1))

        # create group for quantity metadata
        self._qtty_meta = hf.create_group(".quantities")

        # create start iteration
        hiter = hf.create_group("start")
        hiter.attrs["desc"] = "DMFT initial data"
        hf.flush()
        return hiter

    def next_iteration(self, iter_type, iter_no=None):

        if not self.is_writer:
            return

        if iter_no is None:
            self.iter = self.file.create_group(iter_type)
            self.iter.attrs["desc"] = iter_type.capitalize() + " iteration"
        else:
            iter_name = "%s-%03d" % (iter_type.lower(), iter_no + 1)
            self.iter = self.file.create_group(iter_name)
            self.iter.attrs["desc"] = ("%s iteration no. %d"
                                       % (iter_type.capitalize(), iter_no + 1))

            # Set dmft-last etc. to the last iteration other than finish.
            last_link = "%s-last" % iter_type.lower()
            if iter_no >= 1:
                del self.file[last_link]
            current_iter = "/%s-%03d" % (iter_type.lower(), iter_no + 1)
            self.file[last_link] = hdf5.SoftLink(current_iter)

        self.file.flush()
        self.ineq_no = -1
        self.ineq = None

    def _write_meta_return_axes(self, qtty_name):
        try:
            meta_info = meta.QUANTITIES[qtty_name]
        except KeyError:
            warn("No metadata found for quanitity %s." % qtty_name,
                 UserWarning, 2)
            meta_info = {}
        try:
            meta_node = self._qtty_meta.create_group(qtty_name).attrs
        except ValueError:
            pass
        else:
            for key, value in meta_info.iteritems():
                meta_node[key] = value
        try:
            return meta_info["axes"]
        except KeyError:
            return ()

    def write_quantity(self, qtty_name, qtty_data, qtty_shape=None, qtty_slice=None,
                       ineq_force_list=False):
        def write_to_node(node, part, pshape=None):
            if pshape is None:
                try:
                    pshape = part.shape
                except AttributeError:
                    pshape = ()
            pslice = qtty_slice
            if band_dims == 0:
                pass
            elif band_dims == 1:
                # extract the diagonal because that is what should go into HDF5
                part = orbspin.extract_diagonal(part)
                axes_list = [-2, -1] + range(part.ndim - 2)
                part = part.transpose(axes_list)
                pshape = pshape[-2:] + pshape[:-4]
                if pslice is not None:
                    pslice = (slice(None),)*2 + (qtty_slice,)
            elif band_dims == 2:
                axes_list = [-4, -3, -2, -1] + range(part.ndim - 4)
                part = part.transpose(axes_list)
                pshape = pshape[-4:] + pshape[:-4]
                if pslice is not None:
                    pslice = (slice(None),)*4 + (qtty_slice,)
            else:
                warn("Unclear how to deal with > 2 band dims",
                     UserWarning, 3)
            if write_everything:
                node.create_dataset(qtty_name + "/value", data=part)
            else:
                #print "SETTING", qtty_slice, "OF", qtty_name, "SHAPE", pshape
                dataset = node.require_dataset(qtty_name + "/value",
                                               shape=pshape, dtype=part.dtype,
                                               exact=True)
                dataset[pslice] = part

        if not self.is_writer:
            return

        write_everything = qtty_slice is None
        axes = self._write_meta_return_axes(qtty_name)
        band_dims = len({"band", "band1", "band2", "band3", "band4"} & set(axes))
        if axes and axes[0] == "ineq":
            # writing ineq stuff
            axes = axes[1:]
            for iineq, ineq in enumerate(self.ineq_list):
                ineq_node = self.iter.require_group("ineq-%03d" % (iineq+1))
                if ineq_force_list or isinstance(qtty_data, (list, tuple)):
                    # We also expect a shape list here:
                    if qtty_shape is None:
                        ineq_shape = None
                    else:
                        ineq_shape = qtty_shape[iineq]
                    write_to_node(ineq_node, qtty_data[iineq], ineq_shape)
                elif band_dims:
                    # Write quantity that needs to be downfolded
                    if qtty_shape is None:
                        ineq_shape = None
                    else:
                        ineq_shape = ineq.d_downfolded_shape(qtty_shape)
                    write_to_node(ineq_node, ineq.d_downfold(qtty_data),
                                  ineq_shape)
                else:
                    # Write scalar quantity
                    write_to_node(ineq_node, qtty_data)
        else:
            write_to_node(self.iter, qtty_data)
        self.file.flush()
        
    def write_distributed(self, qtty_name, qtty_parts, qtty_shape, 
                          ineq_force_list=False):
        for pslice, pdata in qtty_parts:
            self.write_quantity(qtty_name, pdata, qtty_shape, pslice,
                                ineq_force_list)

    def write_impurity_result(self, iineq, result):
        if not self.is_writer:
            return

        ineq_node = self.iter.require_group("ineq-%03d" % (iineq+1))
        for qtty_name, qtty_value in result.iteritems():
            self._write_meta_return_axes(qtty_name)
            qtty_node = ineq_node.create_group(qtty_name)
            if qtty_value.dtype.fields:
                for fname in qtty_value.dtype.fields:
                    qtty_node.create_dataset(fname, data=qtty_value[fname])
            else:
                qtty_node.create_dataset("value", data=qtty_value)
        self.file.flush()

    
    def write_impurity_component(self, iineq, component):
        if not self.is_writer:
            return

        ineq_node = self.iter.require_group("ineq-%03d" % (iineq+1))
        
        for qtty_name, qtty_value in component.iteritems(): 
           self._write_meta_return_axes(qtty_name.split("/")[0])
           qtty_node = ineq_node.create_group(qtty_name)
           if qtty_value.dtype.fields:
               for fname in qtty_value.dtype.fields:
                   qtty_node.create_dataset(fname, data=qtty_value[fname])
           else:
               qtty_node.create_dataset("value", data=qtty_value)

        self.file.flush()



