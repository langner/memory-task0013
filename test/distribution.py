#!/usr/bin/env python
# @CulgiVersion: 6.0.0
import os 			
import math 			
import glob 			
from math import * 			
import random 			
import re 			
import datetime 			
import CulgiPy 			
from CulgiPy import CulgiError 			
Culgi = CulgiPy.Culgi_GetInstance( ) 			
Palette = Culgi.GetPaletteCmds( ) 			
GraphicsManager = Culgi.GetGraphicsMgrCmds( ) 			
BuildersManager = Culgi.GetBuildersMgrCmds( ) 			
CalculatorsManager = Culgi.GetCalculatorsMgrCmds( ) 			
AnalyzersManager = Culgi.GetAnalyzersMgrCmds( ) 			
PluginsManager = Culgi.GetPluginsMgrCmds( ) 			
UtilitiesManager = Culgi.GetUtilitiesMgrCmds( ) 			
MappersManager = Culgi.GetMappersMgrCmds( )			
SelectorsManager = Culgi.GetSelectorsMgrCmds( )			
MathManager = Culgi.GetMathMgrCmds( )			

def GPERange( initial, end, step, sign):			
	listItems = [] 			
	if type(initial) is int and type(end) is int and type(step) is int :			
		i = 0             
		item = initial 			
		while ( ((sign == "<" or sign == "<=") and item < end  and step > 0) or ( (sign == ">" or sign == ">=") and item > end and step < 0)) :			
			listItems.append( item )			
			i = i + 1             
			item = initial + i* int(step) 			
		if ( ((sign == "<=") and item <= end  and step > 0) or ( (sign == ">=") and item >= end and step < 0)) :			
			listItems.append( item)			
	else :			
		i = 0.0			
		item = initial			
		while ( ((sign == "<" or sign == "<=") and item < end  and step > 0) or ( (sign == ">" or sign == ">=") and item > end and step < 0)) : 			
			listItems.append( item )			
			i = i + 1.0			
			item = initial + i* float(step)			
		if ( ((sign == "<=") and item <= end  and step > 0) or ( (sign == ">=") and item >= end and step < 0)) : 			
			listItems.append( item)			
	return listItems 			

def RemoveFile( fileName ):			
	Culgi.MpiBarrier()			
	mpi = Culgi.GetProcRank()			
	for name in glob.glob(fileName):			
		if os.access(name,os.W_OK):			
			if (mpi==0) :			
				os.remove( name )			
	Culgi.MpiBarrier()			

def RenameFile( oldName, newName):			
	Culgi.MpiBarrier()			
	mpi = Culgi.GetProcRank()			
	if (mpi==0) :			
		os.rename( oldName, newName)			
	Culgi.MpiBarrier()			

def WriteToFile( content, fileID ):			
	mpi = Culgi.GetProcRank()			
	if (mpi==0) :			
		fileID.write( content )			
		fileID.flush()			
	Culgi.MpiBarrier()			

def OpenFile(fileName, access) :			
	if (access=='w') :			
		mpi = Culgi.GetProcRank()			
		if (mpi==0) :			
			return open(fileName, access)			
		else :			
			return 0			
	else :			
		return open(fileName, access)			

def CloseFile (fileID) :			
	if (fileID != 0) :			
		fileID.close()			

def Scan(reString, matchString, *var):			
	result = re.search( reString, matchString)			
	if(result!=None):			
		return result.groups()			
	else :			
		return var			

def IsEOF(file) :			
	line=file.readline()			
	if(line):			
		EOF = 0			
		file.seek(-len(line),1)			
	else:			
		EOF =  1			
	return EOF			

def NewFolder(folder):			
	mpi = Culgi.GetProcRank()			
	if (mpi==0) :			
		os.makedirs(folder)			
	Culgi.MpiBarrier()			


N = 128
D = int(math.sqrt(N * 3.1415*(0.4+0.5)**2 / 0.9))+1
np = Palette.LoadSoftCoreColloid("np.cof")
box = Palette.CreateMesoBox("box", D, D, 1)
box.CopySoftCoreColloids(np, N)
calc = CalculatorsManager.CreateMBFHybridCalculator()
calc.SetSystem(box)

g = GraphicsManager.CreateGraphics()
g.AddViewable(box)

print N, D

for i in range(N):
    npi = box.GetSoftCoreColloidCmds(i)
    npi.SetBeadDisplayRadius("S", 0.0)
    npi.SetDiffusionFactor(0.01)
    npi.SetConstVelocity("Z", 0.0)
    npi.SetConstAngularVelocity("X", 0.0)
    npi.SetConstAngularVelocity("Y", 0.0)
    npi.SetRotationalFactors(0.0, 0.0, 0.05)
    comi = npi.GetCenterOfMass()
    npi.SetCenterOfMass(comi.GetX(), comi.GetY(), 0.5)
    g.AddViewable(npi)

calc.SetTimeStep(0.05)
calc.SetTemperature(1.0)

params = calc.GetParametersCmds()
params_SC = params.GetSCMParametersCmds()
params_SC.SetA("C", "C", 10)
params_SC.SetA("C", "S", 10)
params_SC.SetA("S", "S", 10)

calc.AddGraphics(g)
calc.Run(1000)

D1 = 64
f = 1.0 * D1 / D
box1 = Palette.CreateMesoBox("box1", D1, D1, 1)
box1.CopySoftCoreColloids(np, N)
g1 = GraphicsManager.CreateGraphics()
g1.AddViewable(box1)
for i in range(N):
    comi = box.GetSoftCoreColloidCmds(i).GetCenterOfMass()
    npi = box1.GetSoftCoreColloidCmds(i)
    npi.SetCenterOfMass(f*comi.GetX(), f*comi.GetY(), 0.5)
    npi.SetBeadDisplayRadius("S", 0.0)
    g1.AddViewable(npi)
g1.Render()

Culgi.StartEventLoop( )
