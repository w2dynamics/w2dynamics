# compileall module from cpython 3.1.5 with __future__ imports added
# (included for syntax check test if Python 2 is used)
#
# Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010
# Python Software Foundation; All Rights Reserved

# PYTHON SOFTWARE FOUNDATION LICENSE VERSION 2
# --------------------------------------------

# 1. This LICENSE AGREEMENT is between the Python Software Foundation
# ("PSF"), and the Individual or Organization ("Licensee") accessing and
# otherwise using this software ("Python") in source or binary form and
# its associated documentation.

# 2. Subject to the terms and conditions of this License Agreement, PSF hereby
# grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce,
# analyze, test, perform and/or display publicly, prepare derivative works,
# distribute, and otherwise use Python alone or in any derivative version,
# provided, however, that PSF's License Agreement and PSF's notice of copyright,
# i.e., "Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010
# Python Software Foundation; All Rights Reserved" are retained in Python alone or
# in any derivative version prepared by Licensee.

# 3. In the event Licensee prepares a derivative work that is based on
# or incorporates Python or any part thereof, and wants to make
# the derivative work available to others as provided herein, then
# Licensee hereby agrees to include in any such work a brief summary of
# the changes made to Python.

# 4. PSF is making Python available to Licensee on an "AS IS"
# basis.  PSF MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR
# IMPLIED.  BY WAY OF EXAMPLE, BUT NOT LIMITATION, PSF MAKES NO AND
# DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS
# FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF PYTHON WILL NOT
# INFRINGE ANY THIRD PARTY RIGHTS.

# 5. PSF SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF PYTHON
# FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS
# A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING PYTHON,
# OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.

# 6. This License Agreement will automatically terminate upon a material
# breach of its terms and conditions.

# 7. Nothing in this License Agreement shall be deemed to create any
# relationship of agency, partnership, or joint venture between PSF and
# Licensee.  This License Agreement does not grant permission to use PSF
# trademarks or trade name in a trademark sense to endorse or promote
# products or services of Licensee, or any third party.

# 8. By copying, installing or otherwise using Python, Licensee
# agrees to be bound by the terms and conditions of this License
# Agreement.
#
"""Module/script to "compile" all .py files to .pyc (or .pyo) file.

When called as a script with arguments, this compiles the directories
given as arguments recursively; the -l option prevents it from
recursing into directories.

Without arguments, if compiles all modules on sys.path, without
recursing into subdirectories.  (Even though it should do so for
packages -- for now, you'll have to deal with packages separately.)

See module py_compile for details of the actual byte-compilation.

"""
from __future__ import absolute_import, division, print_function
import locale
import os
import sys
import py_compile
import struct
import imp

__all__ = ["compile_dir","compile_path"]

def compile_dir(dir, maxlevels=10, ddir=None,
                force=0, rx=None, quiet=0):
    """Byte-compile all modules in the given directory tree.

    Arguments (only dir is required):

    dir:       the directory to byte-compile
    maxlevels: maximum recursion level (default 10)
    ddir:      if given, purported directory name (this is the
               directory name that will show up in error messages)
    force:     if 1, force compilation, even if timestamps are up-to-date
    quiet:     if 1, be quiet during compilation

    """
    if not quiet:
        print('Listing', dir, '...')
    try:
        names = os.listdir(dir)
    except os.error:
        print("Can't list", dir)
        names = []
    names.sort()
    success = 1
    for name in names:
        fullname = os.path.join(dir, name)
        if ddir is not None:
            dfile = os.path.join(ddir, name)
        else:
            dfile = None
        if rx is not None:
            mo = rx.search(fullname)
            if mo:
                continue
        if os.path.isfile(fullname):
            head, tail = name[:-3], name[-3:]
            if tail == '.py':
                if not force:
                    try:
                        mtime = int(os.stat(fullname).st_mtime)
                        expect = struct.pack('<4sl', imp.get_magic(), mtime)
                        cfile = fullname + (__debug__ and 'c' or 'o')
                        with open(cfile, 'rb') as chandle:
                            actual = chandle.read(8)
                        if expect == actual:
                            continue
                    except IOError:
                        pass
                if not quiet:
                    print('Compiling', fullname, '...')
                try:
                    ok = py_compile.compile(fullname, None, dfile, True)
                except KeyboardInterrupt:
                    raise KeyboardInterrupt
                except py_compile.PyCompileError as err:
                    if quiet:
                        print('*** Error compiling', fullname, '...')
                    else:
                        print('*** ', end='')
                    # escape non-printable characters in msg
                    msg = err.msg.encode((sys.stdout.encoding
                                          if sys.stdout.encoding is not None
                                          else locale.getpreferredencoding()),
                                         'backslashreplace')
                    msg = msg.decode((sys.stdout.encoding
                                      if sys.stdout.encoding is not None
                                      else locale.getpreferredencoding()))
                    print(msg)
                    success = 0
                except (SyntaxError, UnicodeError, IOError) as e:
                    if quiet:
                        print('*** Error compiling', fullname, '...')
                    else:
                        print('*** ', end='')
                    print(e.__class__.__name__ + ':', e)
                    success = 0
                else:
                    if ok == 0:
                        success = 0
        elif maxlevels > 0 and \
             name != os.curdir and name != os.pardir and \
             os.path.isdir(fullname) and \
             not os.path.islink(fullname):
            if not compile_dir(fullname, maxlevels - 1, dfile, force, rx,
                               quiet):
                success = 0
    return success

def compile_path(skip_curdir=1, maxlevels=0, force=0, quiet=0):
    """Byte-compile all module on sys.path.

    Arguments (all optional):

    skip_curdir: if true, skip current directory (default true)
    maxlevels:   max recursion level (default 0)
    force: as for compile_dir() (default 0)
    quiet: as for compile_dir() (default 0)

    """
    success = 1
    for dir in sys.path:
        if (not dir or dir == os.curdir) and skip_curdir:
            print('Skipping current directory')
        else:
            success = success and compile_dir(dir, maxlevels, None,
                                              force, quiet=quiet)
    return success

def main():
    """Script main program."""
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'lfqd:x:')
    except getopt.error as msg:
        print(msg)
        print("usage: python compileall.py [-l] [-f] [-q] [-d destdir] " \
              "[-x regexp] [directory ...]")
        print("-l: don't recurse down")
        print("-f: force rebuild even if timestamps are up-to-date")
        print("-q: quiet operation")
        print("-d destdir: purported directory name for error messages")
        print("   if no directory arguments, -l sys.path is assumed")
        print("-x regexp: skip files matching the regular expression regexp")
        print("   the regexp is searched for in the full path of the file")
        sys.exit(2)
    maxlevels = 10
    ddir = None
    force = 0
    quiet = 0
    rx = None
    for o, a in opts:
        if o == '-l': maxlevels = 0
        if o == '-d': ddir = a
        if o == '-f': force = 1
        if o == '-q': quiet = 1
        if o == '-x':
            import re
            rx = re.compile(a)
    if ddir:
        if len(args) != 1:
            print("-d destdir require exactly one directory argument")
            sys.exit(2)
    success = 1
    try:
        if args:
            for dir in args:
                if not compile_dir(dir, maxlevels, ddir,
                                   force, rx, quiet):
                    success = 0
        else:
            success = compile_path()
    except KeyboardInterrupt:
        print("\n[interrupt]")
        success = 0
    return success

if __name__ == '__main__':
    exit_status = int(not main())
    sys.exit(exit_status)
