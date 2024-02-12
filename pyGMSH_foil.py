#!/usr/bin/env python3

"""

MIT License

Copyright (c) 2024 Simão Marques

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import gmsh
import numpy as np
import os,sys
import argparse
import csv
import matplotlib.pyplot as plt
from scipy import interpolate

def run_gmsh(surf1,surf2,isSharp,inputs):
    
    # initialize gmsh
    gmsh.initialize()
    
    #number of points for tfi
    npts1       = inputs["n_points_foil_1"]   # trailing edge to max. thickness location
    npts2       = inputs["n_points_foil_2"]   # max. thickness location to leading edge
    npts3       = inputs["n_points_bl"]    # b.l.
    npts4       = inputs["n_points_te"]    # t.e.
    npts5       = inputs["n_points_wake"]    # wake 
    
    #progression rate for b.l. thickness
    blProgRate      = inputs["blProgRate"]
    surfProgRate    = inputs["surfaceProgRate"]
    
    #b.l. thickness region size;
    blThick1 = inputs["blThick_trailingEdge"]
    blThick2 = inputs["blThick_leadingingEdge"]
    
    # add auxiliary geometry points for ff -- index start at 100    
    radius  = inputs["radius"];
    theta   = np.array([45,115,245,-45])
    x       = radius*np.cos(np.deg2rad(theta))
    y       = radius*np.sin(np.deg2rad(theta))
    
    ind     = 100;
    gmsh.model.geo.add_point(0.5, 0, 0, 1.0,ind)
    for i in range(4):
        ind += 1
        gmsh.model.geo.add_point(x[i]+.5,y[i], 0, 0.25*radius,ind)
    
    # set surf1 as upper surface
    if surf1[1,1] < surf2[1,1]:
        tmp     = surf1.copy();
        surf1   = surf2.copy();
        surf2   = tmp.copy();
    
    #trailing edge tangent vector
    v1  = np.array([surf1[-1,0] - surf1[-2,0],surf1[-1,1] - surf1[-2,1]])
    v2  = np.array([surf2[-1,0] - surf2[-2,0],surf2[-1,1] - surf2[-2,1]])
    vte = 0.5*(v1+v2)
    vte = vte/np.linalg.norm(vte)
    
    newInd  = (-surf1[:, 0]).argsort()
    surf1   = surf1[newInd]
    
    lc       = 1e-2    
    ## aerofoil points start at 1000:    
    ind     = ite1 = 1000
    ile     = 1000 + len(surf1)-1
    ymax    = -1.
    ## upper surface
    for x in surf1:
        gmsh.model.geo.add_point(x[0],x[1],0,lc,ind)
        if x[1] > ymax:
            ymax = x[1]
            itop = ind
        ind += 1
    
    #lower surface
    ind     = 1500
    ymin    = 1.
    ist     = 1 # do not repeat leading edge point
    if isSharp:
        iend = len(surf2) - 1   #do not repeat end point
    else:
        iend = len(surf2)   

    for x in surf2[ist:iend,:]:   
        gmsh.model.geo.add_point(x[0],x[1],0,lc,ind)
        if x[1] < ymin:
            ymin = x[1]
            ibot = ind

        ind += 1
    
    # splines for aerofoil surface
    splineList  = []
    for x in range(1000,itop+1): 
        splineList.append(x)
    spline1     = gmsh.model.geo.add_spline(splineList)
    
    splineList  = []
    for x in range(itop,1000+len(surf1)): 
        splineList.append(x)
    spline2     = gmsh.model.geo.add_spline(splineList)
    
    splineList  = [1000+len(surf1)-1]
    for x in range(1500,ibot+1): 
        splineList.append(x)
    spline3     = gmsh.model.geo.add_spline(splineList)
    
    splineList  = []
    
    if isSharp:
        for x in range(ibot,1500+len(surf2)-2): splineList.append(x)
        splineList.append(1000)
        spline4     = gmsh.model.geo.add_spline(splineList)
        #gmsh.model.geo.addCurveLoop([spline1,spline2,spline3,spline4],4000)
        ite2 = ite1
    else:
        for x in range(ibot,1500+len(surf2)-1): splineList.append(x)
        spline4     = gmsh.model.geo.add_spline(splineList)
        ite2        = 1500+len(surf2)-2
        lineTE      = gmsh.model.geo.add_line(ite2,ite1)
        #gmsh.model.geo.addCurveLoop([spline1,spline2,spline3,spline4,lineTE],4000)


    ## offset points for wake/b.l. region
    ## upper surface
    ind = 2000    
    
    #first point(s) is for wake spliting line
    if isSharp:
        gmsh.model.geo.add_point(surf1[0,0]+0.2*vte[0],surf1[0,1]+0.02*vte[1],0,lc,ind)
        wakeEnd = ind;
        ind += 1;
    else:
        wakeEnd = np.array([0,0])
        gmsh.model.geo.add_point(surf1[0,0]+0.2*vte[0],surf1[0,1]+0.02*np.abs(vte[1]),0,lc,ind)
        wakeEnd[0] = ind;
        ind += 1;
        gmsh.model.geo.add_point(surf2[-1,0]+0.2*vte[0],surf1[-1,1]-0.02*np.abs(vte[1]),0,lc,ind)
        wakeEnd[1] = ind;
        ind += 1;

    # trailing edge
    nx0 = surf1[1,1] - surf1[0,1]         #-dy
    ny0 = surf1[0,0] - surf1[1,0]           #dx
    nx  = blThick1*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
    ny  = blThick1*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
    gmsh.model.geo.add_point(surf1[0,0]+nx+0.2*vte[0],surf1[0,1]+ny+0.05*np.abs(vte[1]),0,lc,ind)  # wake
    wakeEndTop = ind
    
    ind += 1
    gmsh.model.geo.add_point(surf1[0,0]+nx,surf1[0,1]+ny,0,lc,ind)      # b.l. region @ t.e.
    iblpt = []
    iblpt.append(ind)
    
    lineWakeTop = gmsh.model.geo.add_line(wakeEndTop,iblpt[0])
    lineBL1     = gmsh.model.geo.add_line(ite1,iblpt[0])
    
    splineList  = []
    splineList.append(ind)
    #mid points
    ind     += 1
    for index in range(1002,itop+1): 
        i   = index - 1000;
        nx0 = surf1[i+1,1] - surf1[i,1]         #-dy
        ny0 = surf1[i,0] - surf1[i+1,0]           #dx
        nx  = blThick1*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
        ny  = blThick1*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
        gmsh.model.geo.add_point(surf1[i,0]+nx,surf1[i,1]+ny,0,lc,ind)        
        splineList.append(ind)
        itopBL  = ind
        ind     += 1
    
    splineBL1 = gmsh.model.geo.add_spline(splineList)
    lineiTop  = gmsh.model.geo.add_line(itop,itopBL)    # midsection
    
    splineList  = []
    splineList.append(itopBL)
    for index in range(itop+1,1000+len(surf1)-1): 
        i       = index - 1000;
        nx0     = surf1[i+1,1] - surf1[i,1]         #-dy
        ny0     = surf1[i,0] - surf1[i+1,0]           #dx
        beta    = (blThick2*(1-surf1[i,1]/surf1[itop- 1000,1]) + blThick1*surf1[i,1]/surf1[itop- 1000,1])
        nx      = beta*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
        ny      = beta*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
        gmsh.model.geo.add_point(surf1[i,0]+nx,surf1[i,1]+ny,0,lc,ind)        
        splineList.append(ind)
        ind     += 1
    
    #leading edge
    nx0 = surf1[ile-ite1,1] - surf1[ile-ite1-1,1]         #-dy
    ny0 = surf1[ile-ite1-1,0] - surf1[ile-ite1,0]           #dx
    
    nx0 = 0.5*(nx0 + surf2[1,1] - surf2[0,1])         #-dy
    ny0 = 0.5*(ny0 + surf2[0,0] - surf2[1,0])       # dx
    
    nx  = blThick2*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
    ny  = blThick2*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
    gmsh.model.geo.add_point(surf1[ile-ite1,0]+nx,surf1[ile-ite1,1]+ny,0,lc,ind)    
    
    splineList.append(ind)
    splineBL2   = gmsh.model.geo.add_spline(splineList)
    lineBL2     = gmsh.model.geo.add_line(ile,ind)      #leading edge

    ##lower surface
    
    # leading edge
    splineList = []
    splineList.append(ind)  
    iblpt.append(ind)  
    ind += 1
    
    ##max thickness
    for index in range(1501,ibot+2): 
        i       = index - 1500;
        nx0     = surf2[i+1,1] - surf2[i,1]         #-dy
        ny0     = surf2[i,0] - surf2[i+1,0]           #dx
        beta    = (blThick2*(1-surf2[i,1]/surf2[ibot- 1500,1]) + blThick1*surf2[i,1]/surf2[ibot- 1500,1])
        nx      = beta*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
        ny      = beta*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
        gmsh.model.geo.add_point(surf2[i,0]+nx,surf2[i,1]+ny,0,lc,ind)        
        splineList.append(ind)
        ibotBL  = ind
        ind     += 1
    
    splineBL3 = gmsh.model.geo.add_spline(splineList)
    lineiBot  = gmsh.model.geo.add_line(ibot,ibotBL)
    
    splineList = []
    splineList.append(ibotBL)
    for index in range(ibot+2,len(surf2)+1500-1): 
        i   = index - 1500;
        nx0 = surf2[i+1,1] - surf2[i,1]         #-dy
        ny0 = surf2[i,0] - surf2[i+1,0]           #dx
        nx  = blThick1*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
        ny  = blThick1*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
        gmsh.model.geo.add_point(surf2[i,0]+nx,surf2[i,1]+ny,0,lc,ind)        
        splineList.append(ind)
        ind     += 1

    ## trailing edge
    nx0 = surf2[-1,1] - surf2[-2,1]         #-dy
    ny0 = surf2[-2,0] - surf2[-1,0]           #dx
    nx  = blThick1*nx0/np.sqrt(nx0*nx0 + ny0*ny0)
    ny  = blThick1*ny0/np.sqrt(nx0*nx0 + ny0*ny0)
    gmsh.model.geo.add_point(surf2[-1,0]+nx,surf2[-1,1]+ny,0,lc,ind)      # b.l. region @ t.e.
    splineList.append(ind)
    iblpt.append(ind)
    
    splineBL4   = gmsh.model.geo.add_spline(splineList)    
    lineBL3     = gmsh.model.geo.add_line(ite2,iblpt[2])

    #lower side wake
    ind += 1
    gmsh.model.geo.add_point(surf2[-1,0] + nx + 0.2,surf2[-1,1] + ny + 0.05*vte[1],0,lc,ind)  # wake
    wakeEndBot  = ind
    lineWakeBot = gmsh.model.geo.add_line(iblpt[-1],wakeEndBot)
    
    # Wake region TFI
    if isSharp:
        #  blup *-----------------------*wakeEndTop
        #       |                       |
        #  ite1 *-----------------------*wakeEnd
        #       |                       |
        # bldwn *-----------------------*wakeEndBot
        lineSplitWake   = gmsh.model.geo.add_line(ite1,wakeEnd)
        lineWake1       = gmsh.model.geo.add_line(wakeEnd,wakeEndTop)
        lineWake2       = gmsh.model.geo.add_line(wakeEnd,wakeEndBot)
        
        gmsh.model.geo.addCurveLoop([lineSplitWake,lineWake1,lineWakeTop,-lineBL1],4001)
        gmsh.model.geo.addCurveLoop([lineSplitWake,lineWake2,-lineWakeBot,-lineBL3],4002)
        
        # set point distributions
        gmsh.model.geo.mesh.setTransfiniteCurve(lineSplitWake, npts5, "Progression", surfProgRate)
        gmsh.model.geo.mesh.setTransfiniteCurve(lineWake2, npts3, "Progression", surfProgRate)
        
        # add aeroil physical group
        gmsh.model.addPhysicalGroup(1, [spline1,spline2,spline3,spline4], 10000,"foil")
        
    else:
        #  blup *-----------------------*wakeEndTop
        #       |                       |
        #  ite1 *-----------------------*wakeEnd[0]
        #       |                       |
        #  ite2 *-----------------------*wakeEnd[1]
        #       |                       |
        # bldwn *-----------------------*wakeEndBot
        lineSplitWake1  = gmsh.model.geo.add_line(ite1,wakeEnd[0])
        lineSplitWake2  = gmsh.model.geo.add_line(ite2,wakeEnd[1])
        lineWake1       = gmsh.model.geo.add_line(wakeEnd[0],wakeEndTop)
        lineWake2       = gmsh.model.geo.add_line(wakeEnd[1],wakeEndBot)
        lineWake3       = gmsh.model.geo.add_line(wakeEnd[0],wakeEnd[1])
        
        gmsh.model.geo.addCurveLoop([lineSplitWake1,lineWake1,lineWakeTop,-lineBL1],4001)
        gmsh.model.geo.addCurveLoop([lineSplitWake1,lineWake3,-lineSplitWake2,lineTE],4002)
        gmsh.model.geo.addCurveLoop([lineSplitWake2,lineWake2,-lineWakeBot,-lineBL3],4003)
        
        # add aeroil physical group
        gmsh.model.addPhysicalGroup(1, [spline1,spline2,spline3,spline4,lineTE], 10000,"foil")
        
        # set point distributions
        gmsh.model.geo.mesh.setTransfiniteCurve(lineSplitWake1, npts5, "Progression", surfProgRate)
        gmsh.model.geo.mesh.setTransfiniteCurve(lineSplitWake2, npts5, "Progression", surfProgRate)
        gmsh.model.geo.mesh.setTransfiniteCurve(lineWake3, npts4, "Bump", .05)
        gmsh.model.geo.mesh.setTransfiniteCurve(lineTE, npts4, "Bump", .005)
        

    # set point distributions
    gmsh.model.geo.mesh.setTransfiniteCurve(lineWake2, npts3, "Progression", surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(lineWake1, npts3, "Progression", surfProgRate)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(lineBL1, npts3, "Progression", blProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(lineBL3, npts3, "Progression", blProgRate)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(lineWakeTop, npts5, "Progression", -surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(lineWakeBot, npts5, "Progression", surfProgRate)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(splineBL1, npts1, "Progression", surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(spline1, npts1, "Progression", surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(lineiTop, npts3, "Progression", blProgRate)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(splineBL2, npts2, "Progression", -surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(spline2, npts2, "Progression", -surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(lineBL2, npts3, "Progression", blProgRate)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(splineBL3, npts2, "Progression", surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(spline3, npts2, "Progression", surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(lineiBot, npts3, "Progression", blProgRate)
    
    gmsh.model.geo.mesh.setTransfiniteCurve(splineBL4, npts1, "Progression", -surfProgRate)
    gmsh.model.geo.mesh.setTransfiniteCurve(spline4, npts1, "Progression", -surfProgRate)

    #======
    
    gmsh.model.geo.synchronize()
    
    
    if isSharp: 
        gmsh.model.geo.addPlaneSurface([4001], 5001)
        gmsh.model.geo.addPlaneSurface([4002], 5002)
        gmsh.model.geo.mesh.setTransfiniteSurface(5001)
        gmsh.model.geo.mesh.setTransfiniteSurface(5002)
    else:
        gmsh.model.geo.addPlaneSurface([4001], 5001)
        gmsh.model.geo.addPlaneSurface([4003], 5003)
        gmsh.model.geo.addPlaneSurface([4002], 5002)
        gmsh.model.geo.mesh.setTransfiniteSurface(5002)
        gmsh.model.geo.mesh.setTransfiniteSurface(5003)

    gmsh.model.geo.mesh.setTransfiniteSurface(5001)
    gmsh.model.geo.mesh.setRecombine(2, 5001)
    gmsh.model.geo.mesh.setRecombine(2, 5002)
    
    if not isSharp: gmsh.model.geo.mesh.setRecombine(2, 5003)

    # foil BL regions
    gmsh.model.geo.addCurveLoop([splineBL1,-lineiTop,-spline1,lineBL1],4004)
    gmsh.model.geo.addCurveLoop([spline2,lineBL2,-splineBL2,-lineiTop],4005)
    gmsh.model.geo.addCurveLoop([spline3,lineiBot,-splineBL3,-lineBL2],4006)
    gmsh.model.geo.addCurveLoop([spline4,lineBL3,-splineBL4,-lineiBot],4007)
    
    gmsh.model.geo.addPlaneSurface([4004], 5004)
    gmsh.model.geo.addPlaneSurface([4005], 5005)
    gmsh.model.geo.addPlaneSurface([4006], 5006)
    gmsh.model.geo.addPlaneSurface([4007], 5007)

    gmsh.model.geo.mesh.setTransfiniteSurface(5004)
    gmsh.model.geo.mesh.setTransfiniteSurface(5005)
    gmsh.model.geo.mesh.setTransfiniteSurface(5006)
    gmsh.model.geo.mesh.setTransfiniteSurface(5007)
    
    gmsh.model.geo.mesh.setRecombine(2, 5004)
    gmsh.model.geo.mesh.setRecombine(2, 5005)
    gmsh.model.geo.mesh.setRecombine(2, 5006)
    gmsh.model.geo.mesh.setRecombine(2, 5007)
    
    gmsh.model.geo.synchronize()
    
    # farfield boundary
    lc2 = 5;
    ff1 = gmsh.model.geo.addCircleArc(101, 100, 102);
    ff2 = gmsh.model.geo.addCircleArc(102, 100, 103);
    ff3 = gmsh.model.geo.addCircleArc(103, 100, 104);
    ff4 = gmsh.model.geo.addCircleArc(104, 100, 101)

    # add aeroil physical group
    gmsh.model.addPhysicalGroup(1, [ff1,ff2,ff3,ff4], 10001,"farfield")
        
    gmsh.model.geo.addCurveLoop([ff1,ff2,ff3,ff4],4008)
    if isSharp: 
        gmsh.model.geo.addCurveLoop([lineWakeTop,splineBL1,splineBL2,splineBL3,splineBL4,lineWakeBot,-lineWake2,lineWake1],4009)
    else:
        gmsh.model.geo.addCurveLoop([lineWakeTop,splineBL1,splineBL2,splineBL3,splineBL4,lineWakeBot,-lineWake2,-lineWake3,lineWake1],4009)
    
    #main mesh surface
    gmsh.model.geo.addPlaneSurface([4008,4009], 5008)

    gmsh.model.geo.synchronize()
    
    
    ## add fields
    if isSharp:
        wakepnt = wakeEnd
    else:
        wakepnt = wakeEnd[0]

    gmsh.model.mesh.field.add("Distance", 1) 
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [100])
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [splineBL1,splineBL2,splineBL3,splineBL4])
        
    gmsh.model.mesh.field.add("Distance", 2) 
    gmsh.model.mesh.field.setNumbers(2, "PointsList", [wakepnt])

    gmsh.model.mesh.field.add("Threshold", 3) 
    gmsh.model.mesh.field.setNumber(3, "InField", 1) 
    gmsh.model.mesh.field.setNumber(3, "SizeMin", 0.5*lc) 
    gmsh.model.mesh.field.setNumber(3, "SizeMax", radius/3.33) 
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.075) 
    gmsh.model.mesh.field.setNumber(3, "DistMax", radius)
    #gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.model.mesh.field.add("Threshold", 4) 
    gmsh.model.mesh.field.setNumber(4, "InField", 2) 
    gmsh.model.mesh.field.setNumber(4, "SizeMin", .5*lc) 
    gmsh.model.mesh.field.setNumber(4, "SizeMax", radius/3.33) 
    gmsh.model.mesh.field.setNumber(4, "DistMin", 0.125) 
    gmsh.model.mesh.field.setNumber(4, "DistMax", radius)

    gmsh.model.mesh.field.add("Threshold", 5)
    gmsh.model.mesh.field.setNumber(5, "InField", 1) 
    gmsh.model.mesh.field.setNumber(5, "SizeMin", 2*lc) 
    gmsh.model.mesh.field.setNumber(5, "SizeMax", radius/3.33) 
    gmsh.model.mesh.field.setNumber(5, "DistMin", 0.75) 
    gmsh.model.mesh.field.setNumber(5, "DistMax", radius)
    
    gmsh.model.mesh.field.add("Min", 6)
    gmsh.model.mesh.field.setNumbers(6, "FieldsList", [3,4,5])
    gmsh.model.mesh.field.setAsBackgroundMesh(6)

    gmsh.model.geo.synchronize()
    # When the element size is fully specified by a mesh size field (as it is in
    # this example), it is thus often desirable to set

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    #note that bug will write this group and boundary elements in SU2 mesh -- bug has been fixed in new release;
    gmsh.model.addPhysicalGroup(2, [5001,5002,5003,5004,5005,5006,5007,5008], 10002)
    
    # Create the relevant Gmsh data structures from Gmsh model.
    gmsh.model.geo.synchronize()
    
    gmsh.option.setNumber("Mesh.Smoothing", 5)
    
    # Generate mesh:
    gmsh.model.mesh.generate()
    
    # Write mesh data:
    gmsh.write("aerofoil.geo_unrolled")
    if isSharp: 
        gmsh.write("mesh_out_s.su2")
    else:
        gmsh.write("mesh_out_b.su2")
    
    # Creates  graphical user interface
    if 'close' not in sys.argv:
        gmsh.fltk.run()
    
    # It finalize the Gmsh API
    gmsh.finalize()

#=======================================================
def get_aerofoil_coordinates(inFile,plotOption = False):
#=======================================================
    """ extract surface mesh points coordinates from coordinate file 
        
            input:  file name
            output: array for side1, points ordered leading edge -> trailing edge
                    array for side2, points ordered leading edge -> trailing edge
                    flag for trailing edge
    """
    
    isSharp = False
    rows    = np.zeros([9999,2])
    
    with open(inFile, 'r', newline='\n') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(), delimiters=', \t')
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        
        ipt     = 0;
        xmin    = 9999.9    #   location of leading points
        xmax    = -9999.9    #   location of leading points
        xstag   = 9999.9    #   location of stagnation point

        for row in reader:
            line = []
            [line.append(x) for x in row if len(x) > 1]
            rows[ipt,0] = float(line[0])
            rows[ipt,1] = float(line[1])

            ## find leading edge
            if rows[ipt,0] <= xmin: 
                xmin    = rows[ipt,0]
                iXmin   = ipt;
            
            ## find trailing edge
            if rows[ipt,0] >= xmax:
                xmax    = rows[ipt,0]
                iXmax   = ipt

            ipt += 1

    rows        = rows[0:ipt,:]
    
    #normalizing aerofoil and setting l.e. @ (0,0)
    rows[:,0]   = rows[:,0] - min(rows[:,0])
    rows[:,1]   = rows[:,1]/max(rows[:,0])
    rows[:,0]   = rows[:,0]/max(rows[:,0])
    
    if iXmin != 0 and np.abs(rows[iXmin,0] - rows[0,0]) < 1e-12:
        data_1   = rows[0:iXmin,:]
        data_2   = rows[iXmin:,:]    
    elif iXmin == 0:
        if iXmax < ipt-1:
            data_1   = rows[0:iXmax+1,:]
            data_2   = rows[iXmax+1:,:]
    elif iXmax == 0:
        data_1   = rows[0:iXmin+1,:]
        data_2   = rows[iXmin:,:]
    else:
        data_1   = rows[0:iXmin+1,:]
        data_2   = rows[iXmin:,:]
    
    splprep1, fp, ier, msg   = interpolate.splprep([data_1[:,0], data_1[:,1]], u=None, k=5, s=1e-9, per=0, full_output=1)
    splprep2, fp, ier, msg   = interpolate.splprep([data_2[:,0], data_2[:,1]], u=None, k=5, s=1e-9, per=0, full_output=1)

    x1,y1 =  interpolate_foil_coordinates(splprep1,du1=1.e-3,du2=1.e-3,resolution=199,N=99,D=20)
    x2,y2 =  interpolate_foil_coordinates(splprep2,du1=1.e-3,du2=1.e-3,resolution=199,N=99,D=20)
    
    if plotOption:
        plt.plot(data_1[:,0],data_1[:,1],'o',label='data1')
        plt.plot(data_2[:,0],data_2[:,1],'*',label='data2')
        
        plt.plot(x1,y1,'r',label='spl1')
        plt.plot(x2,y2,'b',label='spl2')

        plt.legend()
        plt.grid()
        plt.show()
    
    data_1  = data_1[data_1[:, 0].argsort()]
    data_2  = data_2[data_2[:, 0].argsort()]  
    
    #check if t.e. is blunt
    dx = np.abs(data_1[-1,0] - data_2[-1,0])
    dy = np.abs(data_1[-1,1] - data_2[-1,1])
    
    if dx*dx+dy*dy < 1e-7:
        isSharp = True
   
    return data_1,data_2,isSharp

#====================================================================================
def interpolate_foil_coordinates(spline,du1=1.e-3,du2=1.e-3,resolution=99,N=99,D=20):
#====================================================================================

        tck, u  = spline
        xx      = np.linspace(-np.pi+du1,0-du2, resolution)
        uu      = 0.5 + 0.5*np.cos(xx)
        x, y    = interpolate.splev(uu, tck, der=0)

        return x,y


#==========
def main():
#==========

    basePath = os.getcwd()+'/'
    
    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Aerofoil Mesh Generator based on pyGMSH')
    
    # Required positional argument
    parser.add_argument('foilPath', type=str,help='path to coordinate file')
    
    # Optional argument
    parser.add_argument('--plot', type=str,default=False,help='plot aerofoil and splining')

    args = parser.parse_args()
    print("Argument values:")
    
    X1,X2,teFlag = get_aerofoil_coordinates(args.foilPath,args.plot)
    
    parameters = {
        "n_points_foil_1": 129,             # n. of points in first half of aerofoil | te -> midpoint -> le
        "n_points_foil_2": 135,             # n. of points in second half of aerofoil
        "n_points_bl": 29,                  # trailing edge to max. thickness location
        "n_points_te": 29,                  # t.e.
        "n_points_wake": 65,                # n of points in wake
        "radius": 25,                       # distance to farfield
        "blProgRate": 1.5,                  # progression rate for b.l. thickness
        "surfaceProgRate": 1.01,            # chordwise progression rate on surface
        "blThick_trailingEdge": 0.015,      # b.l. thickness region size;
        "blThick_leadingingEdge": 0.006     # b.l. thickness region size;
    }    
    
    run_gmsh(X1,X2,teFlag,parameters)
    
if __name__ == '__main__':
    main()
