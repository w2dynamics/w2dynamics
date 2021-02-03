from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import time
import os
from warnings import warn
import h5py as hdf5
import numpy as np

import w2dyn.auxiliaries as _aux
import w2dyn.auxiliaries.transform as _tf
import w2dyn.dmft.meta as meta
import w2dyn.auxiliaries.statistics as statistics
import w2dyn.dmft.orbspin as orbspin

import sys
ustr = str if sys.version_info >= (3, ) else unicode
h5ustrs = hdf5.special_dtype(vlen=ustr)


def _flat_items(d, prefix=""):
    items = []
    for k in d:
        v = d[k]
        fullk = prefix + k
        if isinstance(v, dict):
            items.extend(_flat_items(v, fullk + "."))
        else:
            items.append((fullk, v))
    return items

class HdfOutput:
    @staticmethod
    def load_old(oldfile, olditer=-1):
        hf = hdf5.File(oldfile, "r")
        file_version = tuple(hf.attrs["outfile-version"])

        if olditer == -1:
            oldit = hf["dmft-last"]
        else:
            oldit = hf["dmft-%03d" % olditer]

        ineqs = [oldit[name] for name in sorted(oldit)
                 if name.startswith("ineq-")]
        if not ineqs:
            raise ValueError("No impurity results found in iteration")

        mu = oldit["mu/value"][()]
        smom_dd = None
        dc_latt = None
        beta = None

        if file_version >= (2, 1):
            siw_dd = [ineq["siw-full/value"][()].transpose(4,0,1,2,3)
                      for ineq in ineqs]
            smom_dd = [ineq["smom-full/value"][()].transpose(4,0,1,2,3)
                       for ineq in ineqs]

            # Fix shape of smom when reading old data files with only
            # one order (i.e. sigma_inf) written
            if any(smom_ineq.shape[0] == 1 for smom_ineq in smom_dd):
                smom_dd = [np.broadcast_to(smom_ineq,
                                           [2] + list(smom_ineq.shape[1:]))
                           for smom_ineq in smom_dd]

            try:
                dc_latt = oldit["dc-latt/value"][()]
            except KeyError:
                pass
        elif file_version >= (2, 0):
            # FIXME: double-counting for old runs
            siw_dd = [orbspin.promote_diagonal(ineq["siw/value"][()].transpose(2,0,1))
                      for ineq in ineqs]
        else:
            # FIXME: double-counting for very old runs
            siw_dd = [orbspin.promote_diagonal(siw_ineq.transpose(2,0,1))
                      for siw_ineq in oldit["siw/value"][()]]

        try:
            beta = hf[".config"].attrs["general.beta"]
        except KeyError:
            pass
        except AttributeError:
            pass

        hf.close()
        return mu, siw_dd, smom_dd, dc_latt, beta

    @staticmethod
    def load_old_kappa(oldfile, olditer=-1, oldtype='dmft', kappa_sign=None):
        hf = hdf5.File(oldfile, "r")
        file_version = tuple(hf.attrs["outfile-version"])

        if olditer == -1:
            oldit = hf[oldtype + "-last"]
            for igrp, grp in ((i, hf[oldtype + "-%03d" % i])
                              for i in range(1, 1000)
                              if (oldtype + "-%03d" % i) in hf):
                if grp == oldit:
                    olditer = igrp
                    break
            else:
                raise ValueError("Malformed file, %s could not be matched to a"
                                 "numbered iteration" % (oldtype + "-last"))
        else:
            oldit = hf[oldtype + "-%03d" % olditer]

        mu = oldit["mu/value"][()]

        def get_nimps(iteration):
            nimps = 0.0

            # FIXME: for anything older than the (2, 2) format I just guessed
            if file_version >= (2, 0):
                ineqs = [iteration[name] for name in sorted(iteration)
                         if name.startswith("ineq-")]
                if not ineqs:
                    raise ValueError("No impurity results found in iteration")

                for occ in (ineq["occ/value"][()] for ineq in ineqs):
                    nimps += np.trace(np.trace(occ, 0, -3, -1), 0, -2, -1)
            else:
                for occ in iteration["occ/value"][()]:
                    nimps += np.trace(np.trace(occ, 0, -3, -1), 0, -2, -1)
            return nimps

        nimps = get_nimps(oldit)

        if olditer > 1:
            preoldit = hf[oldtype + "-%03d" % (olditer - 1)]
            last_mu = preoldit["mu/value"][()]
            last_nimps = get_nimps(preoldit)
            last_mudiff = mu - last_mu
            kappa = (nimps - last_nimps) / (mu - last_mu)
            olditer -= 1
            # try to find kappa with desired sign in earlier
            # iterations if necessary
            while (kappa_sign is not None
                   and olditer > 1
                   and not np.sign(kappa) == np.sign(kappa_sign)):
                oldit = preoldit
                preoldit = hf[oldtype + "-%03d" % (olditer - 1)]
                mu = last_mu
                nimps = last_nimps
                last_mu = preoldit["mu/value"][()]
                last_nimps = get_nimps(preoldit)
                kappa = (nimps - last_nimps) / (mu - last_mu)
                olditer -= 1
            if ((kappa_sign is None or np.sign(kappa) == np.sign(kappa_sign))
                    and last_mu != mu and last_nimps != nimps):
                hf.close()
                return mu, nimps, kappa, last_mudiff

        hf.close()
        return mu, nimps, None, None

    def __init__(self, config, git_revision=None, start_time=None, 
                 ineq_list=None, mpi_comm=None):
        self.quantities_to_write = config["General"]["WriteOnlyQuant"]
        if mpi_comm is not None:
            self.is_writer = mpi_comm.Get_rank() == 0
        else:
            self.is_writer = True

        if self.is_writer:
            # generate prefix-time file name for HDF5 output file
            if start_time is None:
                start_time = time.localtime()
            runstring = time.strftime("%Y-%m-%d-%a-%H-%M-%S", start_time)
            self.filename = "%s-%s.hdf5" % (config["General"]["FileNamePrefix"],
                                            runstring)
            self.filename = mpi_comm.bcast(self.filename, root=0)

            self.file = hdf5.File(self.filename, "w-")
            self.iter = self._ready_file(self.file, config, git_revision,
                                         start_time)
        else:
            self.filename = ""
            self.filename = mpi_comm.bcast(self.filename, root=0)
            self.file = None
            self.iter = None

        # store handles
        self.ineq_list = ineq_list
    
    def _ready_file(self, hf, config, git_revision, start_time):
        # generate HDF5 file and write important metadata
        hf.attrs["outfile-version"] = _aux.OUTPUT_VERSION
        hf.attrs.create("code-version", list(map(ustr, _aux.CODE_VERSION)),
                        dtype=h5ustrs)
        if git_revision is not None:
            hf.attrs["git-revision"] = git_revision
        hf.attrs["run-date"] = time.strftime("%c", start_time)

        # write configureation of current run
        hfgrp = hf.create_group(".config")
        for key, val in _flat_items(config):
            if val is None: continue
            try:
                hfgrp.attrs[key.lower()] = val
            except TypeError:
                hfgrp.attrs.create(key.lower(), val, dtype=h5ustrs)

        # write environment of current run
        hfgrp = hf.create_group(".environment")
        for k in os.environ:
            v = os.environ[k]
            if v is None: continue
            hfgrp.attrs[k] = v

        # write metadata for possible axes objects
        ax = hf.create_group(".axes")
        beta = config["General"]["beta"]
        qcfg = config["QMC"]
        ax.create_dataset("iw", data=_tf.matfreq(beta, "fermi", 2*qcfg["Niw"]))
        ax.create_dataset("pos-iw", data=_tf.matfreq(beta, "fermi", 2*qcfg["Niw"])[qcfg["Niw"]:])
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

    def open_iteration(self, iter_type, iter_no=None):

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

        self.file.flush()
        self.ineq_no = -1
        self.ineq = None

    def close_iteration(self, iter_type, iter_no=None):

        if not self.is_writer:
            return

        # Set dmft-last etc. to the last completed iteration
        last_link = "%s-last" % iter_type.lower()
        if iter_no >= 1:
            del self.file[last_link]
        current_iter = "/%s-%03d" % (iter_type.lower(), iter_no + 1)
        self.file[last_link] = hdf5.SoftLink(current_iter)

        self.file.flush()

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
            for key in meta_info:
                try:
                    meta_node[key] = meta_info[key]
                except TypeError:
                    meta_node.create(key, meta_info[key], dtype=h5ustrs)
        try:
            return meta_info["axes"]
        except KeyError:
            return ()

    def write_quantity(self, qtty_name, qtty_data, qtty_shape=None, qtty_slice=None,
                       ineq_force_list=False):
        if self.quantities_to_write is not None and qtty_name not in self.quantities_to_write:
            return
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
                axes_list = [-2, -1] + list(range(part.ndim - 2))
                part = part.transpose(axes_list)
                pshape = pshape[-2:] + pshape[:-4]
                if pslice is not None:
                    pslice = (slice(None),)*2 + (qtty_slice,)
            elif band_dims == 2:
                axes_list = [-4, -3, -2, -1] + list(range(part.ndim - 4))
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
        if self.is_writer:
            ineq_node = self.iter.require_group("ineq-%03d" % (iineq+1))
        for qtty_name in result:
            if self.quantities_to_write is not None and qtty_name not in self.quantities_to_write:
                continue
            qtty_value = result[qtty_name]
            if self.is_writer:
                self._write_meta_return_axes(qtty_name.split('/')[0])
                qtty_node = ineq_node.require_group(qtty_name)
            try:
                if isinstance(qtty_value, statistics.DistributedSample):
                    mean = qtty_value.mean()
                    stderr = qtty_value.stderr()
                    if self.is_writer:
                        qtty_node.create_dataset('value', data=mean)
                        qtty_node.create_dataset('error', data=stderr)
                elif qtty_value.dtype.fields:
                    for fname in qtty_value.dtype.fields:
                        if self.is_writer:
                            qtty_node.create_dataset(fname, data=qtty_value[fname])
                else:
                    if self.is_writer:
                        qtty_node.create_dataset("value", data=qtty_value)
            except (OSError, RuntimeError, ValueError):
                sys.stderr.write(
                    "\nWARNING: Ignoring field {} for multiple worm components.\n\n".format(qtty_name))
        if self.is_writer:
            self.file.flush()
