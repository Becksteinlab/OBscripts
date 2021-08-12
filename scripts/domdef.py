#!/usr/bin/env python
# Copyright (c) 2010, 2018 Oliver Beckstein
# Released under the MIT license

"""%prog [options] datafile

Read a domain definition data file and write a file defining selection macros.

Definition datafile format:

  # comment and empty lines are ignored
  domain_name   start_resid  end_resid
  ...

  # Simple compound selections: translated to vmd and Charmm
  # (introduce with '@' in first column, no space between name and @)
  #   |  or
  #   &  and
  #   !  not
  #  ()  kept verbatim
  @domain_name  compound_selection
  ...

By convention, the resids are the numbers found in the x-ray structure. Use
--offset to get psf numbering.

Compound selections are very primitive; typically you should *only* use
domain_names that you defined previously and use boolean operators and
parentheses to join them. No residue numbers are translated in compound
selections.
"""

from __future__ import with_statement, print_function, division

from collections import OrderedDict

import numpy



class DomDef(object):
    """Class to handle input/ouput.

    TODO: Rewrite more nicely so that I can derive classes for various
          output codes instead of write_XXX() methods, or at least
          make it more table-driven.
    """
    # ugly & growing... :
    vmdmacro = "atomselect macro %(domname)s {resid %(start_resid)d to %(end_resid)d}"
    vmdcompound = "atomselect macro %(domname)s {%(definition)s}"
    charmmselection = "define %(domname)-8s select resid %(start_resid)3d : %(end_resid)3d  end"
    charmmcompound = "define %(domname)-8s select %(definition)s  end"
    pymolselection = "select %(domname)s, polymer and resi %(start_resid)d-%(end_resid)d"
    pymolcompound = "select %(domname)s, %(definition)s"
    def __init__(self,filename,offset=0):
        self.filename = filename
        self.domains = OrderedDict()
        self.compounds = OrderedDict()
        self.domain_order = []
        self.offset = int(offset)
        self.load(filename)

    def load(self,filename):
        with open(filename,'r') as domdef:
            for line in domdef:
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                fields = line.split()
                if line.startswith('@'):
                    # simple compound command
                    domname, definition = fields[0][1:], " ".join(fields[1:])
                    self.compounds[str(domname)] = definition
                else:
                    # simple selection (macro) definition
                    domname,start_resid,end_resid = fields
                    self.domains[str(domname)] = [self.transform(int(n)) for n in (start_resid, end_resid)]
                    self.domain_order.append(domname)
        # find first and last residue from macros
        resids = numpy.sort(numpy.ravel(self.domains.values()))
        self.first = resids[0]
        self.last = resids[-1]

    def transform(self, resid):
        """Return the resid with offset applied; minimum value is 1."""
        x = resid + self.offset
        if x < 1:
            x = 1
        return x

    def write(self,filename):
        """Write as domdef input file."""
        with open(filename,'w') as domdef:
            domdef.write("""# $Id: domdef.py 1155 2010-05-17 17:15:26Z oliver $
# This file is in a format suitable for scripts/domdef.py
# input = %(filename)r
# offset = %(offset)d

# name      start    end
"""
                         % vars(self))
            for domname,(start_resid,end_resid) in self._ordered_domains():
                domdef.write("%(domname)-11s %(start_resid)5d  %(end_resid)5d\n"% vars())

            # write compounds
            domdef.write("# compound statements\n")
            for domname, definition in self.compounds.items():
                domdef.write("@%(domname)-11s  %(definition)s\n" % vars())

    def write_vmd(self,filename):
        with open(filename,'w') as tcl:
            tcl.write("# $Id" "$\n"
                      """# domain defs from "%(filename)s"\n"""
                      """# offset = %(offset)d\n""" % vars(self))
            for domname,(start_resid,end_resid) in self._ordered_domains():
                tcl.write(self.vmdmacro % vars() + '\n')

            # write compounds
            tcl.write("# compound statements\n")
            for domname, definition in self.compounds.items():
                definition = self._transform('vmd',definition)
                tcl.write(self.vmdcompound % vars() + '\n')

            self._write_vmd_addrep_domains(tcl)

        print("# Load macros with 'source %s' in VMD" % filename)

    def _write_vmd_addrep_domains(self, tcl, **kwargs):
        kwargs.setdefault('colorID', 6)  # silver
        kwargs.setdefault('material', 'AOChalky')
        kwargs.setdefault('representation', 'cartoon')
        tcl.write(("# functions to display representations\n"+
                   "puts {usage: addrep_domains [colorID [material [rep]]]}\n"+
                   "proc addrep_domains {{color %(colorID)d} {material %(material)s} {representation %(representation)s}} {\n"+
                   "  mol color ColorID $color\n"+
                   "  mol material $material\n"+
                   "  mol rep $representation\n") % kwargs)
        tcl.write("  set selections \"%s\"\n" % " ".join(self.domain_order))
        tcl.write("  foreach sel $selections {\n"
                  "    mol selection $sel\n"
                  "    mol addrep top\n"
                  "  }\n")
        tcl.write("}\n")

    def write_bendix(self, filename):
        resids = []
        for domname,(start_resid,end_resid) in self._ordered_domains():
            resids.extend((str(start_resid), str(end_resid)))
        with open(filename, 'w') as bendix:
            bendix.write(" ".join(resids) + "\n")

        print("# Wrote domains to Bendix helix file %r" % filename)

    def write_pymol(self,filename):
        with open(filename,'w') as pml:
            pml.write("# $Id" "$\n"
                     """# domain defs from "%(filename)s"\n"""
                     """# offset = %(offset)d\n""" % vars(self))
            for domname,(start_resid,end_resid) in self._ordered_domains():
                pml.write(self.pymolselection % vars() + '\n')

            # write compounds
            pml.write("# compound statements\n")
            for domname, definition in self.compounds.items():
                definition = self._transform('pymol',definition)
                pml.write(self.pymolcompound % vars() + '\n')

        print("# Load selection with '@%s' in PyMOL" % filename)

    def write_charmm(self,filename):
        with open(filename,'w') as charmm:
            charmm.write("! $Id" "$\n"
                """! domain defs from "%(filename)s"\n"""
                """! offset = %(offset)d\n""" % vars(self))
            for domname,(start_resid,end_resid) in self._ordered_domains():
                charmm.write(self.charmmselection % vars() + '\n')

            # write compounds
            charmm.write("! compound statements\n")
            for domname, definition in self.compounds.items():
                definition = self._transform('charmm',definition)
                charmm.write(self.charmmcompound % vars() + '\n')

        print("""# Add SELECTION definitons to Charmm script with 'stream "%s"'""" % filename)

    def write_xvg(self,filename):
        """Write secondary structure XVG graph."""

        with open(filename,'w') as domdef:
            domdef.write("""# $Id: domdef.py 1155 2010-05-17 17:15:26Z oliver $
# input = %(filename)r
# offset = %(offset)d
"""
                         % vars(self))
            def _write(yval):
                zero = 0
                start_resid = last_start_resid = self.first
                domdef.write("%(start_resid)5d  %(zero)g\n"% vars())
                for domname,(start_resid,end_resid) in self._ordered_domains():
                    if start_resid < last_start_resid:
                        break  # only look at first set of definitions
                    domdef.write("%(start_resid)5d  %(zero)g\n"% vars())
                    domdef.write("%(start_resid)5d  %(yval)g\n"% vars())
                    domdef.write("%(end_resid)5d  %(yval)g\n"% vars())
                    domdef.write("%(end_resid)5d  %(zero)g\n"% vars())
                    last_start_resid = start_resid
                domdef.write("%5d  %g\n" % (self.last, zero))

            _write(0.5)
            domdef.write('&\n')
            _write(-0.5)
        print("# Wrote xvg file with secondary structure graph")

    def _ordered_domains(self):
        for domname in self.domain_order:
            yield domname, self.domains[domname]

    def _transform(self, mode, s):
        """Replace tokens |&! --> or,and,not"""
        # XXX: should be a dict of translation rules
        if mode == "charmm":
            return s.replace('|', ' .or. ').replace('&', ' .and. ').replace('!',' .not. ')
        elif mode == "vmd":
            return s.replace('|', ' or ').replace('&', ' and ').replace('!',' not ')
        elif mode == "pymol":
            return s.replace('|', ' or ').replace('&', ' and ').replace('!',' not ')
        else:
            raise ValueError('mode = %s not known' % mode)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage=__doc__)
    parser.add_option('-f','--file',dest="filename",
                      help="output domain definition datafile FILE [%default]",metavar="FILE")
    parser.add_option('-t','--vmd',dest="vmdfilename",
                      help="VMD tcl output, defining macro selections",metavar="FILE")
    parser.add_option('-b','--bendix',dest="bendixfilename",
                      help="helix definition file for Bendix",metavar="FILE")
    parser.add_option('-p','--pymol',dest="pymolfilename",
                      help="PyMOL output, defining selections",metavar="FILE")
    parser.add_option('-c','--charmm',dest="charmmfilename",
                      help="Charmm output, defining selections",metavar="FILE")
    parser.add_option('-x','--xvg',dest="xvgfilename",
                      help="XVG output of secondary structure blocks",metavar="FILE")
    parser.add_option('-n','--offset',dest='offset',
                      help="add OFFSET to the resids in the domain file [%default]",metavar="OFFSET")

    parser.set_defaults(offset=0, filename='/dev/stdout')

    opts,args = parser.parse_args()

    try:
        ddfile = args[0]
    except IndexError:
        raise ValueError("domain file is required, see --help")

    D = DomDef(ddfile,offset=opts.offset)

    if opts.filename:
        D.write(opts.filename)
    if opts.vmdfilename:
        D.write_vmd(opts.vmdfilename)
    if opts.bendixfilename:
        D.write_bendix(opts.bendixfilename)
    if opts.pymolfilename:
        D.write_pymol(opts.pymolfilename)
    if opts.charmmfilename:
        D.write_charmm(opts.charmmfilename)
    if opts.xvgfilename:
        D.write_xvg(opts.xvgfilename)
