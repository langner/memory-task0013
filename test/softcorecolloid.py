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



np = Palette.LoadSoftCoreColloid("np.cof")
bcp = Palette.CreateGaussianChain("bcp", "A8B8", 10, 10, 10)
box = Palette.CreateMesoBox("box", 10, 10, 10)
box.AddGaussianChain(bcp)
box.CopySoftCoreColloids(np, 1)
col_0 = box.GetSoftCoreColloidCmds(0)
col_0.SetDiffusionFactor(1)
MBFHybridCalculator = CalculatorsManager.CreateMBFHybridCalculator()
MBFHybridCalculator.SetSystem(box)
MBFHybridCalculator.SetSaveInstantResultsOn("SCMBondEnergy,SCMBendingEnergy,SCMTorsionEnergy,SCMNBEnergy,SCMElectrostaticsEnergy,SCMExternalPotentialEnergy,SCMPotentialEnergy,GCDensityInhomogeneity,GCTotalFreeEnergy,GCIdealFreeEnergy,GCContactFreeEnergy,GCCompressibilityFreeEnergy,GCElectrostaticFreeEnergy,CouplingEnergy")
Culgi.StartEventLoop( )
