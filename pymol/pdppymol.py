import subprocess
import tempfile
import os

def pdp(target, option=""):

#make temporary dir and do everything there
    with tempfile.TemporaryDirectory() as dname:
        #turn off zooming when loading: set auto_zoom, off    
        old_auto_zoom=cmd.get("auto_zoom")
        cmd.set("auto_zoom","off")
        # print tmp dir name
        print("Temporary directory =" + dname)
        # make sure you have pdp_2023 in PATH
        # directly giving 'execute' full path below is good alternative
        # For example : execute = "/usr/bin/pdp++"
        execute = "pdp++"
        tmptarget = dname + "/target.cif"
#        pml = dname + "./naive.pml"
        pml = "./naive.pml"
        # save pdb for mican
        pymol.cmd.save(tmptarget, target)

        mican = [execute, tmptarget]
        
        #for op in option.split():
        #    if(op == "-O"):
        #        print("option -O is reserved")
        #        raise CmdException
        #    mican.append(op)
                
        proc=subprocess.run(mican,stdout = subprocess.PIPE)
        print(proc.stdout.decode("utf8")) # print result to pymol console
        
        pymol.cmd.load(pml)
        color_obj(0)
        cmd.util.cbc(target)
        # reset auto_zoom as you had set
        cmd.set("auto_zoom",old_auto_zoom)

pymol.cmd.extend("pdp",pdp)
cmd.auto_arg[0]['pdp'] = cmd.auto_arg[0]['delete']


def color_obj(rainbow=0):
   """

AUTHOR

    Gareth Stockwell

USAGE

    color_obj(rainbow=0)

    This function colours each object currently in the PyMOL heirarchy
    with a different colour.  Colours used are either the 22 named
    colours used by PyMOL (in which case the 23rd object, if it exists,
    gets the same colour as the first), or are the colours of the rainbow

SEE ALSO

    util.color_objs()
    """

   # Process arguments
   rainbow = int(rainbow)

   # Get names of all PyMOL objects
   obj_list = cmd.get_names("objects")

   if rainbow:

      print("\nColouring objects as rainbow\n")

      nobj = len(obj_list)

      # Create colours starting at blue(240) to red(0), using intervals
      # of 240/(nobj-1)
      for j in range(nobj):
         hsv = (240 - j * 240 / (nobj - 1), 1, 1)
         # Convert to RGB
         rgb = hsv_to_rgb(hsv)
         # Define the new colour
         cmd.set_color("col" + str(j), rgb)
         print(obj_list[j], rgb)
         # Colour the object
         cmd.color("col" + str(j), obj_list[j])

   else:

      print("\nColouring objects using PyMOL defined colours\n")

      # List of available colours
      colours = [
         "red",
         "green",
         "blue",
         "yellow",
         "violet",
         "cyan",
         "salmon",
         "lime",
         "pink",
         "slate",
         "magenta",
         "orange",
         "marine",
         "olive",
         "purple",
         "teal",
         "forest",
         "firebrick",
         "chocolate",
         "wheat",
         "white",
         "grey",
      ]
      ncolours = len(colours)

      # Loop over objects
      i = 0
      for obj in obj_list:
         print("  ", obj, colours[i])
         cmd.color(colours[i], obj)
         i = i + 1
         if i == ncolours:
            i = 0

def tcdock_set_chain_by_cluster():
   for i in range(1, 100):
      print(i)
      cmd.alter("*F_%i_*" % i, "chain=" + ROSETTA_CHAINS[(i - 1) % len(ROSETTA_CHAINS)])
      cmd.alter("*F_%i_*" % i, "chain=" + ROSETTA_CHAINS[(100 - i) % len(ROSETTA_CHAINS)])

# HSV to RGB routine taken from Robert L. Campbell's color_b.py script
#   See http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/
# Original algorithm from: http://www.cs.rit.edu/~ncs/color/t_convert.html

def hsv_to_rgb(hsv):

   h = float(hsv[0])
   s = float(hsv[1])
   v = float(hsv[2])

   if s == 0:
      # achromatic (grey)
      r = g = b = v

   else:
      # sector 0 to 5
      h = h / 60.0
      i = int(h)
      f = h - i  # factorial part of h
      # print h,i,f
      p = v * (1 - s)
      q = v * (1 - s * f)
      t = v * (1 - s * (1 - f))

      if i == 0:
         (r, g, b) = (v, t, p)
      elif i == 1:
         (r, g, b) = (q, v, p)
      elif i == 2:
         (r, g, b) = (p, v, t)
      elif i == 3:
         (r, g, b) = (p, q, v)
      elif i == 4:
         (r, g, b) = (t, p, v)
      elif i == 5:
         (r, g, b) = (v, p, q)
      else:
         (r, g, b) = (v, v, v)
         print("error, i not equal 1-5")

   return [r, g, b]

# Add color_obj to the PyMOL command list
cmd.extend("color_obj", color_obj)
