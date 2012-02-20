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



N = 32
D = 16
box = Palette.CreateMesoBox("box", D, D, 1)
bcp = Palette.CreateGaussianChain("bcp", "A32B20", D, D, 1)
bcp.CreateHomogeneousRelativeDensity(1)
bcp_A = bcp.GetBeadTypeRelativeDensityFieldCmds("A")
np = Palette.LoadSoftCoreColloid("np.cof")
box.AddGaussianChain(bcp)
box.CopySoftCoreColloids(np, N)
Graphics = GraphicsManager.CreateGraphics()
Graphics.AddViewable(box)
Graphics.AddViewable(bcp_A)
for i in GPERange( 0, N, 1, "<" ) :
	col_i = box.GetSoftCoreColloidCmds(i)
	col_i_GeometricCenter = col_i.GetGeometricCenter()
	X = col_i_GeometricCenter.GetX()
	Y = col_i_GeometricCenter.GetY()
	col_i.SetGeometricCenter(X, Y, 0.5)
	col_i.SetConstAngularVelocity(0, 0, 0)
#	col_i.SetConstVelocity(0, 0, 0)
	col_i.SetConstVelocity("Z", 0)
	col_i.SetDiffusionFactor(0.01)
	Graphics.AddViewable(col_i)

MBFHybridCalculator = CalculatorsManager.CreateMBFHybridCalculator()
MBFHybridCalculator.SetSystem(box)
MBFHybridCalculator.SetFieldDiffusionMethod("EntropyFieldDynamics")
MBFHybridCalculator.SetSaveInstantResultsOn("SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy")
MBFHybridCalculator.SetTemperature(0.1)
MBFHybridCalculator_Parameters = MBFHybridCalculator.GetParametersCmds()
MBFHybridCalculator_Parameters.SetBeadFieldCoupling("C", "A", 1)
MBFHybridCalculator_Parameters.SetBeadFieldCoupling("C", "B", 2)
MBFHybridCalculator_Parameters.SetBeadFieldCoupling("S", "A", 1)
MBFHybridCalculator_Parameters.SetBeadFieldCoupling("S", "B", 2)
MBFHybridCalculator_Parameters_SCMParameters = MBFHybridCalculator_Parameters.GetSCMParametersCmds()
MBFHybridCalculator_Parameters_SCMParameters.SetA("S", "S", 25.0)
MBFHybridCalculator_Parameters_SCMParameters.SetA("S", "C", 25.0)
MBFHybridCalculator_Parameters_SCMParameters.SetA("C", "C", 25.0)
MBFHybridCalculator_Parameters_GCParameters = MBFHybridCalculator_Parameters.GetGCParametersCmds()
MBFHybridCalculator_Parameters_GCParameters.SetChi("A", "B", 0.5)
MBFHybridCalculator_Parameters_GCParameters.SetExpansionParameter("A", 0.1)
MBFHybridCalculator_Parameters_GCParameters.SetExpansionParameter("B", 0.1)
ArchivePlugin = PluginsManager.CreateArchivePlugin()
MBFHybridCalculator.AddPlugin(ArchivePlugin)
MBFHybridCalculator.SetGraphicsUpdateFrequency(5)
Graphics.Render()
MBFHybridCalculator.AddGraphics(Graphics)
MBFHybridCalculator.Run(100000000)
Culgi.StartEventLoop( )
