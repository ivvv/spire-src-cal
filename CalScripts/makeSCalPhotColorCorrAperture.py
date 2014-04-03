#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2014 Herschel Science Ground Segment Consortium
#
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
# 
#  Herschel-SPIRE Colour Correction factors 
# 
#  This routine calculates colour correction tables to convert from the standard
#  pipeline flux densities (which are quoted for a nu*F_nu=const. spectrum) into 
#  flux densities for a range of source spectra. These produce monochromatic flux
#  densities at the SPIRE reference wavelengths of 250, 350 and 500 micron. The
#  source spectra include power law spectra, and modified black body spectra with
#  a range of emissivities. Metadata values are produced for the flux conversion
#  parameters applied to the pipeline (the "K4P" and "K4E" parameters).
#
#  Input:
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles and beam model parameters from SPIRE calibration tree
#    Flux conversion parameters from SPIRE calibration tree
#    Name and version of output file
#
#  Output:
#    SCalPhotColorCorrAperture product
#
#===============================================================================
# 
#  Herschel-SPIRE Aperture Corrections
# 
#  This routine calculates the colour corrections for a range of source spectra,
#  and for both point sources and extended sources
# 
#  Input:
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles & metadata from SPIRE calibration tree
#    Example beam map (for header, metadata etc.)
#    Name and version of output file
# 
#  Output:
#    SCalPhotColorCorrAperture_noBG product
#    SCalPhotColorCorrAperture_incBG product
# 
#  Calculations:
#   1. Reads in beam profiles and RSRF from calibration tree, including the
#      values of gamma and effective frequencies
#   2. Scales the monochromatic beam with frequency by exponent gamma and
#      calculates the monochromatic beam area
#        Beam_mono(r,nu) = Beam_full(r*(nu/nu_eff)^-gamma)
#           (meaning that the FWHM varies as (nu/nu_eff)^gamma
#   3. Synthesises the effective beam profile for given source
#      spectrum (currently only for alpha)
#   4. Performs aperture photometry on the source, both including and excluding
#      the background annulus
#   5. Divides the aperture photometry result by the total beam area to give
#      the aperture correction
#
#  There are some functions defined (common between calibration products)
#   * spireMonoBeam: Calculate monochromatic beam profile & area at a given freq
#   * spireMonoAreas: Calculate monochromatic beam areas at range of frequencies
#   * spireEffBeam: Calculate the effective beam profile, area and beam map
#
#===============================================================================
# $Id: makeSCalPhotColorCorrAperture.py,v 1.5 2014/02/27 10:01:13 epoleham Exp $
# 
#  Edition History
#   E. Polehampton   22-10-2013  - First version adapted from Andreas' script - SPCAL-83
#   E. Polehampton   31-10-2013  - Correction to columns read from csv file.
#   E. Polehampton   21-01-2014  - Update descriptions and values - SPCAL-94
#  Chris North   - 18/Feb/2014 - 1.0: First version (full calculation)
#                                     Used for spire_cal_12_1
#                                     SPCAL-109
#  Chris North   - 26/Feb/2014        Corrected calls to spireEffBeam and spireEffBeamMap
#                                     Removed analytical method section
#===============================================================================
import os
scriptVersionString = "makeSCalPhotColorCorrAperture.py $Revision: 1.5 $"

#-------------------------------------------------------------------------------
#===============================================================================
#=====                              DIRECTORIES                            =====
#===============================================================================
#-------------------------------------------------------------------------------

directory = "..//..//..//..//..//..//data//spire//cal//SCal"
dataDir = "//disks//winchester2//calibration_data//"
#LOCAL VERSION
#directory = Configuration.getProperty('var.hcss.workdir')
#dataDir = Configuration.getProperty('var.hcss.workdir')

# if inputCalDirTree is True, then calibration products are read from a
#  calibration directory tree
#  e.g. RSRF is read from <directory>/Phot/SCalPhotRsrf/SCalPhotRsrf.fits
# If inputCalDirTree is False, then required calibration products are read from
#  a calibration file (either a local pool or a jarfile)
inputCalDirTree=True
if not inputCalDirTree:
	#read calibration tree from pool
	cal=spireCal(pool='spire_cal_12_1')
	## alternatively, read from jarFile
	#cal=spireCal(pool=os.path.join(dataDir,'spire_cal_12_1_test.jar'))

# if outputCalDirTree is True, then calibration products are written to a
#  calibration directory tree
#  e.g. <directory>/Phot/SCalPhotColorCorrAperture/SCalPhotColorCorrAperture_noBG_<version>.fits
# If outputCalDirTree is False, then calibration products are written to
#  dataDir e.g. <dataDir>/SCalPhotColorCorrAperture_noBG_<version>.fits
outputCalDirTree=True

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           INPUT PARAMETERS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

# Colour correction table version
version = "3"

# set format version and date format
formatVersion = "1.0"
df  = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")

# set verbose to print more information during processing
verbose = True

#-------------------------------------------------------------------------------
# Input parameters for aperture correction
apPhotRad={"PSW":22.,"PMW":30.,"PLW":45.}
apPhotBGRad={'in':60.,'out':90.}
# range of alphas to compute colour corrections for
alphaK=[-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
#alphaK=[-4,-1,2]

# range of beta and temp to calculate colour corrections for
#betaK=[0.,0.5,1.,1.25,1.5,1.75,2.,2.5,3.]
#tempK=range(3,300)

# Example beam map to read in (to copy header info etc.)
beamMapName="0x5000241aL_PSW_pmcorr_1arcsec_norm_beam.fits"

#-------------------------------------------------------------------------------
#===============================================================================
#=====                         READ FILTER PROFILES                        =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Loading physical constants
from herschel.share.unit import *
h = Constant.H_PLANCK.value
k = Constant.K_BOLTZMANN.value
c = Constant.SPEED_OF_LIGHT.value

# Frequency raster of common frequency grid
deltaNu = 0.1e9		# 0.1 GHz
nuMin   = 300.e9
nuMax   = 1800.e9
nNu     = FIX((nuMax-nuMin)/deltaNu)
freq    = Double1d(range(nNu)) * deltaNu + nuMin

#-------------------------------------------------------------------------------
# Load SPIRE filter functions

# SPIRE band names
spireBands=["PSW","PMW","PLW"]

# Three SPIRE filter reference frequencies for PSW, PMW, PLW respectively
spireRefWl = {"PSW":250.*1e-6, "PMW":350.*1.e-6, "PLW":500.*1.e-6}
spireRefFreq = {}
for band in spireBands:
	spireRefFreq[band] = c/spireRefWl[band]

if inputCalDirTree:
	# SPIRE Photometer RSRF calibration product from cal directory tree
	rsrfVersion = "2"
	rsrf = fitsReader("%s//Phot//SCalPhotRsrf//SCalPhotRsrf_v%s.fits"%(directory, rsrfVersion))
	if verbose:
		print 'Reading RSRF version %s from calibration directory tree'%(rsrfVersion)
else:
	#read RSRF from calibration tree
	rsrf=cal.getPhot().getProduct('Rsrf')
	rsrfVersion=rsrf.getVersion()
	if verbose:
		print 'Reading RSRF version %s from calibration %s'%(rsrfVersion,cal.getVersion())

# Photometer RSRF
spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz
#indexes of freq in rsrf
ixR = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# spire RSRF only
spireFiltOnly={}

#interpolate to freq array
for band in spireBands:
	#create Rsrf object
	interpRsrf = LinearInterpolator(spireFreq, rsrf.getRsrf(band))
	#make arrays for final object
	spireFiltOnly[band] = Double1d(nNu)
	#interpolate Rsrf to freq array
	spireFiltOnly[band][ixR] = interpRsrf(freq[ixR])

#-------------------------------------------------------------------------------
# Load SPIRE Beam Color Corrections
kBeamVersion = "3"
kBeam=fitsReader("%s//Phot//SCalPhotColorCorrBeam//SCalPhotColorCorrBeam_v%s.fits"%(directory, kBeamVersion))

#-------------------------------------------------------------------------------
# Load SPIRE Beam profiles
beamProfsVersion = "3"
beamProfs = fitsReader("%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits"%(directory, beamProfsVersion))
spireEffFreq = {"PSW":beamProfs.meta['freqEffPsw'].double*1.e9,\
	"PMW":beamProfs.meta['freqEffPmw'].double*1.e9,\
	"PLW":beamProfs.meta['freqEffPlw'].double*1.e9}

#-------------------------------------------------------------------------------
arcsec2Sr = (Math.PI/(60.*60.*180))**2
spireAreaEffFreq = {"PSW":beamProfs.meta['beamNeptunePswArc'].double,\
	"PMW":beamProfs.meta['beamNeptunePmwArc'].double,\
	"PLW":beamProfs.meta['beamNeptunePlwArc'].double}
alphaNep={"PSW":beamProfs.meta['alphaNeptunePsw'].double,\
	"PMW":beamProfs.meta['alphaNeptunePsw'].double,\
	"PLW":beamProfs.meta['alphaNeptunePsw'].double}
# Exponent of powerlaw describing FWHM dependence on frequency
# FWHM ~ freq**gamma
gamma = beamProfs.meta['gamma'].double

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           DEFINE FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

# Set up functions to calculate beam profile, effective frequency and effective area
# Function list:
#   * spireMonoBeam: Calculate monochromatic beam profile and area at a given frequency
#   * spireMonoAreas: Calculate monochromatic beam areas over a range of frequencies
#   * hpXcalKcorr: Calculate K-correction parameters for given spectrum & source type

#-------------------------------------------------------------------------------
# Calculate monochromatic beam profile and area at a given frequency
def spireMonoBeam(freqx,beamRad,beamProfs,beamConst,effFreq,gamma,array):
	"""
	========================================================================
	Implements the full beam model to generate the monochromatic beam profile
	and corresponding monochromatic beam solid angle at a given frequency.

	Inputs:
	  freqx:     (float) frequency [Hz] for which monochromatic beam
	               profile and area should be calculated
	  beamRad:   (array float) radius vector from the input beam profiles
	  beamProfs: (dataset) PhotRadialCorrBeam object
	               [used to retrieve core beam profile]
	  beamConst: (array float) constant part of beam profile for "array"
                       [passed to prevent repeated calls]
	  effFreq:   (float) effective frequency [Hz] of array
	  gamma:     (float) Exponent of powerlaw describing FWHM dependence
	               on frequency
	  array:     (string) spire array ('Psw'|'Pmw'|'Plw')

	Outputs:     (list of objects)
	  [0]:       (float) Beam area [arcsec^2] at frequency freqx
	  [1]:       (array float) Monochromatic beam profile at frequency freqx

	Calculation:
          Scales the core beam profile width as (freqx/effFreq)^gamma.
	  Queries the calibration file to generate new core beam profile.
	  Uses constant beam profile where it is larger than core profile.
	  Integrates over radius to caluclate beam area.

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator
	  
	2013/12/19  C. North  initial version

	"""

	#calculate the "scaled" radius, as nu^-gamma
	radNew=beamRad*(freqx/effFreq)**-gamma
	maxRad=max(beamRad)
	nRad=len(beamRad)
	#ensure it doesn't go out of range
	radNew[radNew.where(radNew > maxRad)]=maxRad
	#get the corresponding beam profiles
	beamNew=Double1d(nRad)
	for r in range(nRad):
		beamNew[r]=beamProfs.getCoreCorrection(radNew[r],array)
	#apply the "constant" beam where appropriate
	#beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data
	isConst=beamNew.where(beamNew < beamConst)
	beamNew[isConst]=beamConst[isConst]

	#integrate to get solid angle (in arcsec^2)
	
	beamInterp=LinearInterpolator(beamRad,beamNew * 2. * Math.PI * beamRad)
	integrator=TrapezoidalIntegrator(0,maxRad)
	beamMonoArea=integrator.integrate(beamInterp)

	return(beamMonoArea,beamNew)

#-------------------------------------------------------------------------------
# Calculate monochromatic beam areas over a range of frequencies
def spireMonoAreas(freq,beamProfs,effFreq,gamma,array,freqFact=100):

	"""
	========================================================================
	Generates array of monochromatic beam areas over frequency range by
	calculating over a sparser array and interpolating

	Inputs:
	  freq:      (array float) frequency vector [Hz] for which monochromatic
	               beams areas should be calculated
	  beamProfs: (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreq:   (float) effective frequency [Hz] of array
	  gamma:     (float) Exponent of powerlaw describing FWHM dependence
	               on frequency
	  array:     (string) spire array ('Psw'|'Pmw'|'Plw')
	  freqFact:  (int) Factor by which to reduce size of freq.
	               OPTIONAL. Default=100.

	Outputs:     
	             (array float) Monochromatic Beam area [sr] at frequencies
	                corresponding to freq

	Calculation:
          Geneates sparse frequency array of full range
	  Uses spireMonoBeam to calculate monochromatic beam area at sparse freqs
	  Interpolates to full frequency grid

	Dependencies:
	  spireMonoBeam
	  herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
	  
	2013/12/19  C. North  initial version

	"""

	#set up a sparser range of frequencies (otherwise it takes too long)
	nNu=len(freq)
	nNuArea=nNu/freqFact + 1
	#array of indices of full frequency array to use
	iNuArea=Int1d(range(nNuArea))*freqFact
	iNuArea[-1]=nNu-1

	#set up arrays
	beamMonoFreqSparse=Double1d(nNuArea)
	beamMonoAreaSparse=Double1d(nNuArea)

	#get beam radius array from calibration table
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	#get constant beam profile from calibration table
	beamConst=beamProfs.getConstantCorrectionTable().getColumn(array).data

	# calculate at sparse frequencies
	for fx in range(nNuArea):
		#get corresponding index in full frequency array
		f=iNuArea[fx]
		#populate frequency array
		beamMonoFreqSparse[fx]=freq[f]
		#populate beam area array
		beamMonoAreaSparse[fx]=spireMonoBeam(freq[f],beamRad,beamProfs,beamConst,effFreq,gamma,array)[0]

	# interpolate to full frequency array and convert to Sr
	beamInterp=CubicSplineInterpolator(beamMonoFreqSparse,beamMonoAreaSparse)
	beamMonoArea=beamInterp(freq)*arcsec2Sr #in sr
	
	return(beamMonoArea)

#-------------------------------------------------------------------------------
# Calculate the effective beam profile, area and beam map for SPIRE
def spireEffBeam(freq, transm, beamProfs, effFreq, gamma, array,
  BB=False,temp=20.0,beta=1.8,alpha=-1.0,verbose=False):

	"""
	========================================================================
	Computes an effective beam profile for a given source spectrum
	***
	N.B. Computing the beam area this way integrates over frequency *then* radius,
	  while spireEffArea integrates over radius then frequency.
	  The method used here produces areas which area lower by ~0.1%
	***

	Inputs:
	  freq:       (array float) frequency vector [Hz] for which monochromatic
	                beams areas should be calculated
	  transm:     (array float) relative spectral response (RSRF) corresponding to freq
	                Note that this should *not* include the aperture efficiency
	  beamProfs:  (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreq:    (float) effective frequency [Hz] of array
	  gamma:      (float) Exponent of powerlaw describing FWHM dependence
	                on frequency
	  array:      (string) spire array ('Psw'|'Pmw'|'Plw')
	  BB:         (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	                OPTIONAL. Default=False
	  temp:       (float) Dust/sky temperature (if BB=True)
	                OPTIONAL. Default=20.0
	  beta:       (float) Dust/sky spectral index (if BB=True)
	                OPTIONAL. Default=1.8
	  alpha:      (float) Exponent of power-law sky background model (if BB=False)
	                OPTIONAL. Default=-1

	Outputs:
	              (array float) radialised beam profile

	Calculation:
	  Calculates the source spectrum (either modifies black body or power law)
	  Loops over beam radius, and calculates scaled radii for full frequency range
	  For that radius, gets the values of profile at those scaled radii
	  Integrates profile values over frequency, weighted by RSRF and source spectrum

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2014/01/16  C. North  initial version

	"""

	#
	# Calculate sky background model
	#
	if BB == 1:
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha

	#integrate transm*fSky over frequency for nomalisation
	integrator=TrapezoidalIntegrator(min(freq),max(freq))
	denomInterp=CubicSplineInterpolator(freq,transm*fSky)
	denomInteg=integrator.integrate(denomInterp)

	#get beam radius list from calibration table
	beamRad=beamProfs.getCoreCorrectionTable().getColumn('radius').data
	#get core beam profile from calibration table
	beamCore=beamProfs.getCoreCorrectionTable().getColumn(array).data
	#create interpolation object
	beamCoreInt=CubicSplineInterpolator(beamRad,beamCore)

	#make array for new beam
	nRad=len(beamRad)
	maxRad=max(beamRad)
	effBeam=Float1d(beamRad)
	#loop over radius
	for r in range(nRad):
		#calculate the "scaled" radius for range of frequencies
		radFreq=beamRad[r]*(freq/effFreq)**-gamma
		#ensure it doesn't fo beyong maximum radius
		radFreq[radFreq.where(radFreq > maxRad)]=maxRad
		#compute value beam profile at each scaled radius
		beamCoreFreq=beamCoreInt(radFreq)
		#apply constant beam profile value where appropriate
		beamConstRad=beamProfs.getConstantCorrection(beamRad[r],array)
		isConst=beamCoreFreq.where(beamCoreFreq < beamConstRad)
		beamCoreFreq[isConst]=beamConstRad

		#integrate beamCoreFreq*transm*fSky over frequency
		numInterp=CubicSplineInterpolator(freq,beamCoreFreq*transm*fSky)
		numInteg = integrator.integrate(numInterp)

		#write value into table
		effBeam[r]=numInteg/denomInteg	

	beamInterp = CubicSplineInterpolator(beamRad,effBeam * 2.*Math.PI*beamRad)
	effBeamAreaInt=TrapezoidalIntegrator(0,maxRad)
	effBeamArea=effBeamAreaInt.integrate(beamInterp)
	effBeamDict = {'area':effBeamArea,'profile':effBeam}

	return(effBeamDict)

#-------------------------------------------------------------------------------
# Calculate the effective beam profile, area and beam map for SPIRE
def spireEffBeamMap(beamRad,effBeam,beamRadMap,verbose=False):

	"""
	========================================================================
	Computes an effective beam profile for a given source spectrum
	***
	N.B. Computing the beam area this way integrates over frequency *then* radius,
	  while spireEffArea integrates over radius then frequency.
	  The method used here produces areas which area lower by ~0.1%
	***

	Inputs:
	  beamRad:    (array float) radius array for beam profile [arcsec]
	  effBeam:    (array float) radial beam profile
	  beamRadMap: (Simple Image) 2D image containing radius at each point (in arcsec)

	Outputs:
	              (Simple Image) 2D image containing synthesised beam profile
	                at each point

	Calculation:
	  Creates 2D version of beam profile, based on an input map containing
           the radius from the beam centre at each point

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.interp.CubicSplineInterpolator

	2014/01/16  C. North  initial version
	2014/02/20  C. North  Now only makes 2d beam map
	                      (beam profile is in spireEffBeam)

	"""

	#create beam image by copying beamRad
	effBeamMap=beamRadMap.copy()
	#remove data in image
	effBeamMap.setUnit('Jy/beam')
	effBeamMap['image'].data[:,:]=0

	effBeamInterp=LinearInterpolator(beamRad,effBeam)
	nxMap=int(effBeamMap['image'].data[:,0].size)
	nyMap=int(effBeamMap['image'].data[0,:].size)
	maxRad=max(beamRad)
	if verbose:
		print 'making beam map'
	for x in range(nxMap):
		for y in range(nyMap):
			if beamRadMap['image'].data[x,y] <= maxRad-1:
				effBeamMap['image'].data[x,y]= \
				    effBeamInterp(beamRadMap['image'].data[x,y])

	return(effBeamMap)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           END OF FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#===============================================================================
#=====               CALCULATE APERTURE CORRECTIONS (ANALYTICAL)           =====
#===============================================================================
#-------------------------------------------------------------------------------
##
## set up product for no background
## define the start and end dates for the product
## starting at beginning of PFM1
#startDate = df.parse("2005.02.22/00:00:00/GMT")
#endDate   = df.parse("2020.01.01/00:00:00/GMT")
#
#apCorrNoBG = herschel.spire.ia.dataset.PhotColorCorrAperture() 
#apCorrNoBG.setDescription('SPIRE aperture correction product')
#apCorrNoBG.meta["creator"].value = scriptVersionString
#apCorrNoBG.meta["modelName"].value = "FM"
#apCorrNoBG.meta["creationDate"].value = FineTime(java.util.Date())
#apCorrNoBG.meta["startDate"].value = FineTime(startDate)
#apCorrNoBG.meta["endDate"].value = FineTime(endDate)
#apCorrNoBG.meta["author"]  = herschel.ia.dataset.StringParameter(value="Chris North", description="Author of the data")
##apCorrNoBG.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="%s"%inputFileName, description="Origin of the data")
#apCorrNoBG.meta["dependency"].value = "apertureCorrectionType"
#apCorrNoBG.meta["apertureCorrectionType"] = StringParameter(value="noBG", description="")
#apCorrNoBG.setVersion(version)
#apCorrNoBG.setFormatVersion(formatVersion)
#
##copy product for including background
#apCorrIncBG = apCorrNoBG.copy()
#apCorrNoBG.meta["apertureCorrectionType"] = StringParameter(value="incBG", description="")
#
## Create tables for alpha arrays for aperture corrections
#apCorrNoBG['alpha']=TableDataset(description='Aperture Correction Analytical results without background (Spectral Index)')
#apCorrIncBG['alpha']=TableDataset(description='Aperture Correction Analytical results including background (Spectral Index)')
#apCorrNoBG['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
#apCorrIncBG['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
#
## Add columns for aperture corrections
#for band in spireBands:
#	apCorrNoBG['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
#	apCorrIncBG['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
#
#
#beamRad = beamProfs['core']['radius'].data
#sizeInterp = CubicSplineInterpolator(beamRad,2.*Math.PI*beamRad)
#integBG = TrapezoidalIntegrator(apPhotBGRad['in'],apPhotBGRad['out'])
#integTot = TrapezoidalIntegrator(0,max(beamRad))
#sizeBG = integBG.integrate(sizeInterp)
#
#for band in spireBands:
#	integAp = TrapezoidalIntegrator(0,apPhotRad[band])
#	sizeAp = integAp.integrate(sizeInterp)
#	apCorrNoBG.meta['sizeAp_%s'%band]=DoubleParameter(sizeAp)
#	apCorrNoBG.meta['sizeBG_%s'%band]=DoubleParameter(sizeBG)
#	apCorrIncBG.meta['sizeAp_%s'%band]=DoubleParameter(sizeAp)
#	apCorrIncBG.meta['sizeBG_%s'%band]=DoubleParameter(sizeBG)
#
#print '\nCalculating aperture corrections corrections over alpha (analytical method)...'
#for band in spireBands:
#	for a in range(len(alphaK)):
#		effBeam_x=spireEffBeam(freq,spireFiltOnly[band],beamProfs,spireEffFreq[band],\
#		  gamma,band,BB=False,alpha=alphaK[a],verbose=verbose)
#		print '%s alpha=%.1f'%(band,alphaK[a])
#		# interpolate beam profile
#		beamInterp = CubicSplineInterpolator(beamRad,\
#				effBeam_x['profile'] * 2.*Math.PI*beamRad)
#		# calculate beam areas in aperture and background annulus
#		omegaAp = integAp.integrate(beamInterp)
#		omegaBG = integBG.integrate(beamInterp)
#		# calculate aperture corrections analytically
#		apCorrNoBG['alpha'][band].data[a] = \
#			effBeam_x['area']/omegaAp
#		apCorrIncBG['alpha'][band].data[a] =  \
#			effBeam_x['area']/(omegaAp - omegaBG*sizeAp/sizeBG)
#
## set FITS filenames
#if outputCalDirTree:
#	apCorrNoBGFits = java.io.File(r"%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_Analytical_noBG_v%s.fits"%(directory, version))
#	apCorrIncBGFits = java.io.File(r"%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_Analytical_incBG_v%s.fits"%(directory, version))
#else:
#	apCorrNoBGFits = java.io.File(r"%s//SCalPhotColorCorrAperture_noBG_v%s.fits"%(dataDir, version))
#	apCorrIncBGFits = java.io.File(r"%s//SCalPhotColorCorrAperture_incBG_v%s.fits"%(dataDir, version))
#apCorrNoBG.meta['fileName'] = herschel.ia.dataset.StringParameter(value=apCorrNoBGFits.name,\
#  description="Name of file when exported")
#apCorrIncBG.meta['fileName'] = herschel.ia.dataset.StringParameter(value=apCorrIncBGFits.name,\
#  description="Name of file when exported")
#
##-----------------------------------------------------------------------
## write to FITS files
#fitsWriter = FitsArchive()
#fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
#fitsWriter.save(apCorrNoBGFits.toString(), apCorrNoBG)
#
#fitsWriter = FitsArchive()
#fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
#fitsWriter.save(apCorrIncBGFits.toString(), apCorrIncBG)

#-------------------------------------------------------------------------------
#===============================================================================
#=====              CALCULATE APERTURE CORRECTIONS (FULL METHOD)           =====
#===============================================================================
#-------------------------------------------------------------------------------

#read in example map from file
try:
	beamIn = fitsReader(file = os.path.join(dataDir,beamMapName))
except:
	#download if not available
	import urllib
	urllib.urlretrieve ("https://nhscsci.ipac.caltech.edu/spire/data/beam_profiles/"+beamMapName,\
	    os.path.join(dataDir,beamMapName))
	beamIn = fitsReader(file = os.path.join(dataDir,beamMapName))

beamRad = beamProfs['core']['radius'].data

# make map of radius (speeds up processing later)
print 'Making beam Radius map'
beamRadMap=SimpleImage()
beamRadMap['image']=beamIn['image']
nxMap=beamIn['image'].data[:,0].size
nyMap=beamIn['image'].data[0,:].size
bcenter=[int(nxMap/2.),int(nyMap/2.)]
for x in range(nxMap):
	for y in range(nyMap):
		beamRadMap['image'].data[x,y]= \
			Math.sqrt((x-bcenter[0])**2 + (y-bcenter[1])**2)

# set up product for no background
# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

apCorrFullNoBG = herschel.spire.ia.dataset.PhotColorCorrAperture() 
apCorrFullNoBG.setDescription('SPIRE aperture correction product')
apCorrFullNoBG.meta["creator"].value = scriptVersionString
apCorrFullNoBG.meta["modelName"].value = "FM"
apCorrFullNoBG.meta["creationDate"].value = FineTime(java.util.Date())
apCorrFullNoBG.meta["startDate"].value = FineTime(startDate)
apCorrFullNoBG.meta["endDate"].value = FineTime(endDate)
apCorrFullNoBG.meta["author"]  = herschel.ia.dataset.StringParameter(value="Chris North", description="Author of the data")
#apCorrFullNoBG.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="%s"%inputFileName, description="Origin of the data")
apCorrFullNoBG.meta["dependency"].value = "apertureCorrectionType"
apCorrFullNoBG.meta["apertureCorrectionType"] = StringParameter(value="noBG", description="")
apCorrFullNoBG.setVersion(version)
apCorrFullNoBG.setFormatVersion(formatVersion)

#copy product for including background
apCorrFullIncBG = apCorrFullNoBG.copy()
apCorrFullIncBG.meta["apertureCorrectionType"] = StringParameter(value="incBG", description="")

#-----------------------------------------------------------------------
# Compute aperture corrections for range of alpha
# Create tables for alpha arrays for aperture corrections
apCorrFullNoBG['alpha']=TableDataset(description='Aperture Correction without background (Spectral Index)')
apCorrFullIncBG['alpha']=TableDataset(description='Aperture Correction including background (Spectral Index)')
apCorrFullNoBG['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
apCorrFullIncBG['alpha'].addColumn('alpha',Column(Float1d(alphaK)))

# Add columns for aperture corrections
for band in spireBands:
	apCorrFullNoBG['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
	apCorrFullIncBG['alpha'].addColumn(band,Column(Float1d(len(alphaK))))

#-----------------------------------------------------------------------
print '\nCalculating aperture corrections corrections over alpha...'
for band in spireBands:
	for a in range(len(alphaK)):
		print '%s alpha=%.1f'%(band,alphaK[a])
		#calculate beam areas
		#effBeamPSW contains (beam area, beam profile, beam map)
		effBeam_x=spireEffBeam(freq,spireFiltOnly[band],beamProfs,spireEffFreq[band],\
		  gamma,band,BB=False,alpha=alphaK[a],verbose=verbose)
		effBeam_x['map']=spireEffBeamMap(beamRad,effBeam_x['profile'],beamRadMap,verbose=True)
		#perform aperture photometry
		apPhot_x = annularSkyAperturePhotometry(image=effBeam_x['map'], \
		  fractional=0, centerX=bcenter[0], centerY=bcenter[1], \
		  radiusArcsec=apPhotRad[band], \
		  innerArcsec=apPhotBGRad['in'], outerArcsec=apPhotBGRad['out'])
		#get result of aperture correction procedure
		apPhotIncBG_x = apPhot_x.getTargetTotal()
		apPhotNoBG_x = apPhot_x.getTargetPlusSkyTotal()
		#compute correction to provide actual beam area
		apCorrFullIncBG['alpha'][band].data[a]=effBeam_x['area']/apPhotIncBG_x
		apCorrFullNoBG['alpha'][band].data[a]=effBeam_x['area']/apPhotNoBG_x

# set FITS filenames
if outputCalDirTree:
	apCorrFullNoBGFits = java.io.File(r"%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_noBG_v%s.fits"%(directory, version))
	apCorrFullIncBGFits = java.io.File(r"%s//Phot//SCalPhotColorCorrAperture//SCalPhotColorCorrAperture_incBG_v%s.fits"%(directory, version))
else:
	apCorrFullNoBGFits = java.io.File(r"%s//SCalPhotColorCorrAperture_Full_noBG_v%s.fits"%(dataDir, version))
	apCorrFullIncBGFits = java.io.File(r"%s//SCalPhotColorCorrAperture_Full_incBG_v%s.fits"%(dataDir, version))

apCorrFullNoBG.meta['fileName'] = herschel.ia.dataset.StringParameter(value=apCorrFullNoBGFits.name,\
  description="Name of file when exported")
apCorrFullIncBG.meta['fileName'] = herschel.ia.dataset.StringParameter(value=apCorrFullIncBGFits.name,\
  description="Name of file when exported")

#-----------------------------------------------------------------------
# write to FITS files
fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(apCorrFullNoBGFits.toString(), apCorrFullNoBG)

fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(apCorrFullIncBGFits.toString(), apCorrFullIncBG)
