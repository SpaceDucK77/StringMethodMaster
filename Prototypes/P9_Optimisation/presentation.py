"""
>>> import VIS_col
>>> b = VIS_col.VIS_collection("parameters.smg")
>>> b.parse_parameters()
>>> make_index(b.start, b.index)
>>> sCV = a.startVIS.get_CVs()
get_angle:
                       :-) GROMACS - gmx gangle, 2020.1 (-:

                            GROMACS is written by:
     Emile Apol      Rossen Apostolov      Paul Bauer     Herman J.C. Berendsen
    Par Bjelkmar      Christian Blau   Viacheslav Bolnykh     Kevin Boyd
 Aldert van Buuren   Rudi van Drunen     Anton Feenstra       Alan Gray
  Gerrit Groenhof     Anca Hamuraru    Vincent Hindriksen  M. Eric Irrgang
  Aleksei Iupinov   Christoph Junghans     Joe Jordan     Dimitrios Karkoulis
    Peter Kasson        Jiri Kraus      Carsten Kutzner      Per Larsson
  Justin A. Lemkul    Viveca Lindahl    Magnus Lundborg     Erik Marklund
    Pascal Merz     Pieter Meulenhoff    Teemu Murtola       Szilard Pall
    Sander Pronk      Roland Schulz      Michael Shirts    Alexey Shvetsov
   Alfons Sijbers     Peter Tieleman      Jon Vincent      Teemu Virolainen
 Christian Wennberg    Maarten Wolf      Artem Zhmurov
                           and the project leaders:
        Mark Abraham, Berk Hess, Erik Lindahl, and David van der Spoel

Copyright (c) 1991-2000, University of Groningen, The Netherlands.
Copyright (c) 2001-2019, The GROMACS development team at
Uppsala University, Stockholm University and
the Royal Institute of Technology, Sweden.
check out http://www.gromacs.org for more information.

GROMACS is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 2.1
of the License, or (at your option) any later version.

GROMACS:      gmx gangle, version 2020.1
Executable:   /usr/local/gromacs/bin/gmx
Data prefix:  /usr/local/gromacs
Working dir:  /Users/markopetrovic/Box Sync/Egenstudier/Åk5/P3/Ex-Jobb/StringMethodMaster/Prototypes/Test_P9
Command line:
  gmx gangle -g1 dihedral -group1 dihedrals -s '/Users/markopetrovic/Box Sync/Egenstudier/Åk5/P3/Ex-Jobb/StringMethodMaster/Prototypes/Test_P9/mdrun_2_302333861_0/confout.gro' -n index.ndx -oall temp.xvg


-------------------------------------------------------
Program:     gmx gangle, version 2020.1
Source file: src/gromacs/commandline/cmdlineparser.cpp (line 275)
Function:    void gmx::CommandLineParser::parse(int *, char **)

Error in user input:
Invalid command-line options
  In command-line option -n
    File 'index.ndx' does not exist or is not accessible.
    The file could not be opened.
      Reason: No such file or directory
      (call to fopen() returned error code 2)

For more information and tips for troubleshooting, please check the GROMACS
website at http://www.gromacs.org/Documentation/Errors
-------------------------------------------------------

>>> sCV = a.startVIS.get_CVs()
>>> co = np.array(a.CV_targets)
>>> eCV = a.endVIS.get_CVs()
>>> eCV
[179.918, 70.748, -69.388, -179.468]
>>> co
array([[-179.2324,  -66.8773,   59.0429,  179.5609],
       [-179.3268,  -51.5856,   44.7728,  179.6688],
       [-179.4212,  -36.2939,   30.5027,  179.7767],
       [-179.5156,  -21.0022,   16.2326,  179.8846],
       [-179.61  ,   -5.7105,    1.9625,  179.9925],
       [-179.7044,    9.5812,  -12.3076, -179.8996],
       [-179.7988,   24.8729,  -26.5777, -179.7917],
       [-179.8932,   40.1646,  -40.8478, -179.6838],
       [-179.9876,   55.4563,  -55.1179, -179.5759]])
>>> co=np.concatenate((np.array(sCV[1:3]),co[:,1:3]),axis = 0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<__array_function__ internals>", line 6, in concatenate
ValueError: all the input arrays must have same number of dimensions, but the array at index 0 has 1 dimension(s) and the array at index 1 has 2 dimension(s)
>>> co=np.concatenate((np.array([sCV[1:3]]),co[:,1:3]),axis = 0)
>>> co
array([[-76.379 ,  75.299 ],
       [-66.8773,  59.0429],
       [-51.5856,  44.7728],
       [-36.2939,  30.5027],
       [-21.0022,  16.2326],
       [ -5.7105,   1.9625],
       [  9.5812, -12.3076],
       [ 24.8729, -26.5777],
       [ 40.1646, -40.8478],
       [ 55.4563, -55.1179]])
>>> co=np.concatenate((co[:,1:3], np.array([eCV[1:3]])),axis = 0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<__array_function__ internals>", line 6, in concatenate
ValueError: all the input array dimensions for the concatenation axis must match exactly, but along dimension 1, the array at index 0 has size 1 and the array at index 1 has size 2
>>> co=np.concatenate((co[:,1:3], np.array(eCV[1:3])),axis = 0)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<__array_function__ internals>", line 6, in concatenate
ValueError: all the input arrays must have same number of dimensions, but the array at index 0 has 2 dimension(s) and the array at index 1 has 1 dimension(s)
>>> co=np.concatenate((co, np.array([eCV[1:3]])),axis = 0)
>>> co
array([[-76.379 ,  75.299 ],
       [-66.8773,  59.0429],
       [-51.5856,  44.7728],
       [-36.2939,  30.5027],
       [-21.0022,  16.2326],
       [ -5.7105,   1.9625],
       [  9.5812, -12.3076],
       [ 24.8729, -26.5777],
       [ 40.1646, -40.8478],
       [ 55.4563, -55.1179],
       [ 70.748 , -69.388 ]])
>>> import matplotlib.pyplot as plt
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> plt.plot(co[0,0],co[0,1],"ro")
[<matplotlib.lines.Line2D object at 0x1271fdf10>]
>>> plt.plot(co[-1,0],co[-1,1],"ro")
[<matplotlib.lines.Line2D object at 0x1272085d0>]
>>> plt.plot(co[:,0],co[:,1],"-x")
[<matplotlib.lines.Line2D object at 0x1271f8c10>]
>>> plt.savefig("linear.png")
>>> plt.close()
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> a.strings
[<VIS_col.VIS_string object at 0x123547890>, <VIS_col.VIS_string object at 0x12357eb10>, <VIS_col.VIS_string object at 0x1235b1690>, <VIS_col.VIS_string object at 0x1235e6150>, <VIS_col.VIS_string object at 0x123623bd0>, <VIS_col.VIS_string object at 0x12365b690>, <VIS_col.VIS_string object at 0x12368f150>, <VIS_col.VIS_string object at 0x1236bdbd0>, <VIS_col.VIS_string object at 0x1236f2690>, <VIS_col.VIS_string object at 0x123729150>, <VIS_col.VIS_string object at 0x123756bd0>, <VIS_col.VIS_string object at 0x12378b690>, <VIS_col.VIS_string object at 0x123802150>, <VIS_col.VIS_string object at 0x123831bd0>, <VIS_col.VIS_string object at 0x123864690>, <VIS_col.VIS_string object at 0x12389b150>]
>>> phie
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'phie' is not defined
>>> phie = np.array([sCV[CV1index], eCV[CV1index]])
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'CV1index' is not defined
>>> phie = np.array([sCV[0], eCV[0]])
>>> psie = np.array([sCV[1], eCV[1]])
>>> plt.plot(phie,psie,"ro-")
[<matplotlib.lines.Line2D object at 0x1276bcc50>]
>>> plt.show()
>>> phie = np.array([sCV[1], eCV[1]])
>>> phie = np.array([sCV[2], eCV[2]])
>>> plt.plot(phie,psie,"ro-")
[<matplotlib.lines.Line2D object at 0x121bd8a10>]
>>> plt.show()
>>> plt.plot(phie,psie,"ro-")
[<matplotlib.lines.Line2D object at 0x128fd24d0>]
>>> plt.show()
>>> win = plt.figure()
>>> win
<Figure size 640x480 with 0 Axes>
>>> a.strings.plot_CVs_2D(win,(1,2))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
AttributeError: 'list' object has no attribute 'plot_CVs_2D'
>>> a.strings[2].plot_CVs_2D(win,(1,2))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/Users/markopetrovic/Box Sync/Egenstudier/Åk5/P3/Ex-Jobb/StringMethodMaster/Prototypes/Test_P9/VIS_col.py", line 450, in plot_CVs_2D
    p = plotwindow.plot(self.start_CVs[1:-1, CV_index[0]],
AttributeError: 'Figure' object has no attribute 'plot'
>>> a.strings[2].plot_CVs_2D(plt,(1,2))
>>> plt.savefig("s3orig.png")
>>> plt.xlabel("phi")
Text(0.5, 23.52222222222222, 'phi')
>>> plt.ylabel("psi")
Text(35.472222222222214, 0.5, 'psi')
>>> plt.title("String starting points")
Text(0.5, 1.0, 'String starting points')
>>> plt.savefig("s3orig.png")
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> plt.xlabel("phi")
Text(0.5, 0, 'phi')
>>> plt.ylabel("psi")
Text(0, 0.5, 'psi')
>>> plt.title("String drifted")
Text(0.5, 1.0, 'String drifted')
>>> s = a.strings[2]
>>> plt.plot(s.drift_CVs[1:-1, 1],s.drift[1:-1,2], "bo")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
AttributeError: 'VIS_string' object has no attribute 'drift'
>>> plt.plot(s.drift_CVs[1:-1, 1],s.drift_CVs[1:-1,2], "bo")
[<matplotlib.lines.Line2D object at 0x121dbc350>]
>>> plt.savefig("s3drifted.png")
>>> plt.savefig("s3orig_drifted.png")
>>> plt.plot(s.start_CVs[:, 1],s.start_CVs[:,2], "b")
[<matplotlib.lines.Line2D object at 0x121dc6b90>]
>>> plt.savefig("s3orig_drifted.png")
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> plt.xlabel("phi")
Text(0.5, 0, 'phi')
>>> plt.ylabel("psi")
Text(0, 0.5, 'psi')
>>> plt.title("String drifted")
Text(0.5, 1.0, 'String drifted')
>>> plt.title("String spline")
Text(0.5, 1.0, 'String spline')
>>> plt.plot(s.drift_CVs[1:-1, 1],s.drift_CVs[1:-1,2], "bo")
[<matplotlib.lines.Line2D object at 0x1224a8350>]
>>> plt.savefig("s3spline_drifted.png")
>>> np
<module 'numpy' from '/Users/markopetrovic/Box Sync/Egenstudier/Åk5/P3/Ex-Jobb/StringMethodMaster/masters_venv/lib/python3.7/site-packages/numpy/__init__.py'>
>>> t = np.linspace(0,len(s.SO)+1, 200)
>>> splines = s.spline_data[0]
>>> mins = s.spline_data[1]
>>> deltas = s.spline_data[2]
>>> xdim = 1
>>> ydim = 2
>>> x = splines[xdim](t) * deltas[xdim] + mins[xdim]
>>> y = splines[ydim](t) * deltas[ydim] + mins[ydim]
>>> plt.plot(x,y,"b")
[<matplotlib.lines.Line2D object at 0x121e05810>]
>>> plt.savefig("s3spline_drifted.png")
>>> plt.savefig("s3spline_drifted_targets.png")
>>> plt.plot(s.new_CVs[1:-1, 1],s.new_CVs[1:-1,2], "bv")
[<matplotlib.lines.Line2D object at 0x1224e2fd0>]
>>> plt.savefig("s3spline_drifted_targets.png")
>>> plt.close()
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> plt.close("all")
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> plt.xlabel("phi")
Text(0.5, 0, 'phi')
>>> plt.ylabel("psi")
Text(0, 0.5, 'psi')
>>> plt.title("Steered MD")
Text(0.5, 1.0, 'Steered MD')
>>> plt.savefig("s3steered.png")
>>> plt.plot(s.new_CVs[1:-1, 1],s.new_CVs[1:-1,2], "bv")
[<matplotlib.lines.Line2D object at 0x123703b90>]
>>> plt.savefig("s3steered.png")
>>> plt.plot(s.start_CVs[:, 1],s.start_CVs[:,2], "b")
[<matplotlib.lines.Line2D object at 0x1236f2590>]
>>> plt.savefig("s3steered.png")
>>> plt.close()
>>> plt.figure()
<Figure size 640x480 with 0 Axes>
>>> plt.xlabel("phi")
Text(0.5, 0, 'phi')
>>> plt.ylabel("psi")
Text(0, 0.5, 'psi')
>>> plt.title("Updates")
Text(0.5, 1.0, 'Updates')
>>> iter = 1
>>> plt.title("Iteration: " + str(iter))
Text(0.5, 1.0, 'Iteration: 1')
>>> plt.plot(s.start_CVs[:, 1],s.start_CVs[:,2], "b")
[<matplotlib.lines.Line2D object at 0x1224d7790>]
>>> plt.savefig("sanim_"+str(iter)+".png")
>>> def myplot(s,iter):
...  plt.figure()
...  plt.xlabel("phi")
...  plt.ylabel("psi")
...  plt.title("Iteration: " + str(iter))
...  plt.plot(s.start_CVs[:, 1],s.start_CVs[:,2], "b")
...  plt.savefig("sanim_"+str(iter)+".png")
...
>>> myplot(a.strings[5], 1)
>>> plt.close("all")
>>> for i, ns in enumerate(a.strings):
...  myplot(ns, i)
...  plt.close()
...
>>> for i, ns in enumerate(a.strings):
...  myplot(ns, i+1)
...  plt.close()
...
>>>
"""
