#!/usr/bin/env python3
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from configobj import ConfigObj
from configobj import Section
import optparse
import os
import shutil

class Cfg_mod(ConfigObj):
   """class to modify neseted config objects"""
   def __init__(self,cfg):
      self.cfg=cfg
   
   def read_from_section(self,section_list,key):
      """reads key from nested cfg"""
      tmp=self.cfg 
      for section in section_list:
         tmp=tmp[section]
      
      if key in tmp:
         return tmp[key]
      else:
         raise KeyError("Key not found in Section_List")
      
   def write_to_section(self,section_list,key,value):
      """writes key,value pair to nested cfg"""
      tmp = {key:value}
      
      for section in section_list[::-1]:
        tmp = {section:tmp}
      
      self.cfg.merge(tmp)
      
if __name__ == "__main__":
   """This program serves as a simple converter, from old config files to new nested config files."""

   
   parser = optparse.OptionParser(version="1.0",
                usage="""%prog [--help] [--cfg <cfg infile>] [--cspc <cfg_spec file>] [--o <cfg outfile>]""",
                description="""Converts old Parameters File to new Parameters File using nesting
                defined in Configspec File""")
   parser.add_option("--cfg", metavar="File", dest="cfg_file",
                     help="Input File, default is Parameters.in")
   parser.add_option("-o", metavar="File", dest="outfile",
                     help="Output File, default is to replace Input File and save a copy to Parameters_old.in")
   parser.add_option("--cspc", metavar="File", dest="cspc_file",
                     help="ConfigSpec to evaluate nesting, default is ./auxiliaries/configspec")
   
   options, filenames = parser.parse_args()
   
   if options.cfg_file is not None:
      print("Reading Config File " + options.cfg_file)
      cfg_infile = ConfigObj(options.cfg_file,file_error=True)
      if cfg_infile is None:
         raise IOError("Cannot find Config File ") + options.cfg_file
   else:
      print("Reading default Config File Parameters.in")
      cfg_infile = ConfigObj("Parameters.in",file_error=True)
      if cfg_infile is None:
         raise IOError("Cannot find default Config File Parameters.in")
   
   if options.cspc_file is not None:
      print("Reading ConfigSpec File " + options.cspc_file)
      config_spec = ConfigObj(options.cspc_file,file_error=True,list_values=False,_inspec=True)
      if config_spec is None:
         raise IOError("Cannot find ConfigSpec File ") + options.cspc_file
   else:
      print("Reading default ConfigSpec File ./auxiliaries/configspec")
      #default path to configpsec
      path = os.path.dirname(__file__) + "/auxiliaries/configspec"
      config_spec = ConfigObj(path,file_error=True,list_values=False,_inspec=True)
      if config_spec is None:
         raise IOError("Cannot find default ConfigSpec File ./auxiliaries/configspec")
   
   
   #load sections from config_spec to evaluate nesting up to depth 3, only keys rel vals set to zero
   section_list={}
   for section in config_spec:
      if isinstance(config_spec[section],Section):
         section_list[section] = {}
         for sub_section in config_spec[section]:
            if isinstance(config_spec[section][sub_section],Section):
               section_list[section][sub_section] = {}
               for key in config_spec[section][sub_section]:
                  section_list[section][sub_section][key] = 0
            else:
               section_list[section][sub_section] = 0
      else:
         section_list[section] = 0
  
   #NAt is necessary to determine how many ineqs need to be created
   try: 
      cfg_infile["NAt"]
   except Exception:
      raise ValueError("Please set NAt in Input File")
   
   cfg_outfile = ConfigObj()
   cfg_out = Cfg_mod(cfg_outfile)
   
   for key in cfg_infile:
      if key in section_list["Atoms"]["__many__"]:
         for i in range(int(cfg_infile["NAt"])):
            cfg_out.write_to_section(("Atoms",str(i+1)),key,cfg_infile[key])
      elif key in section_list["QMC"]:
         cfg_out.write_to_section(("QMC",),key,cfg_infile[key])
      elif key in section_list["General"]:
         cfg_out.write_to_section(("General",),key,cfg_infile[key])
      else:
         raise KeyError("Key not in configspec")
   
   if options.outfile is not None:
      print("Creating Output File " + options.outfile)
      cfg_outfile.filename=options.outfile
      if cfg_outfile is None:
         raise IOError("Cannot create Output File ") + options.outfile
   else:
      if options.cfg_file is not None:
         shutil.copy2("./"+options.cfg_file,"./Parameters_old.in")
         cfg_outfile.filename=options.cfg_file
         if cfg_outfile is None:
            raise IOError("Cannot create Output File ") + options.outfile
      else:
         print("Writing to standard Output File")
         shutil.copy2("./Parameters.in","./Parameters_old.in")
         cfg_outfile.filename="Parameters.in"
         if cfg_outfile is None:
            raise IOError("Cannot create Output File ") + options.outfile
   
   #writing data to file
   cfg_out.cfg.write()
   
