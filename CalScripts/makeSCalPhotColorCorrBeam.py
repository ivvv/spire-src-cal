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
#    RSRF from SPIRE calibration tree
#    Beam profiles and beam model parameters from SPIRE calibration tree
#    Name and version of output file
#
#  Output:
#    SCalPhotColorCorrBeam product
#
#===============================================================================
# 
#  Herschel-SPIRE Beam Colour Corrections
# 
#  This routine calculates the beam correction for a range of source spectra.
# 
#  Input:
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles & metadata from SPIRE calibration tree
#    Neptune spectral index and measured area
#    Array of alphas and (beta,temperatures) from SPIRE calibration tree
#    Name and version of output file
# 
#  Output:
#    SCalPhotColorCorrBeam product
# 
#  Calculations:
#   1. Reads in beam profiles and RSRF from calibration tree, including the
#      values of gamma, effective beam area and the pipeline beam areas
#   2. Scales the monochromatic beam with frequency by exponent gamma and
#      calculates the monochromatic beam area
#        Beam_mono(r,nu) = Beam_full(r*(nu/nu_eff)^-gamma)
#           (meaning that the FWHM varies as (nu/nu_eff)^gamma
#        Omega_mono(nu) = int{Beam_mono(r,nu) * 2*pi*r dr}
#   3. By integrating over the RSRF and source spectrum, it calculates the
#      effective beam area for the pipeline spectrum
#        Omega_eff(alpha) = int{Omega_mono(nu) * rsrf(nu) * nu^alpha dnu} / 
#                                 int{rsrf(nu) * nu^alpha dnu}
#   4. Computes the beam correction Kbeam as
#       Kbeam(alpha) = Omega_eff(-1) / Omega_eff(alpha)
#      And similarly for greybody spectra for a range of betas and temperature
#
#  There are some functions defined (common between calibration products)
#   * spireMonoBeam: Calculate monochromatic beam profile & area at a given freq
#   * spireMonoAreas: Calculate monochromatic beam areas at range of frequencies
#   * spireEffArea: Calculate the effective beam area for a given spectrum
#
#===============================================================================
# $Id: makeSCalPhotColorCorrBeam.py,v 1.6 2014/02/26 15:40:46 epoleham Exp $
# 
#  Edition History
#   E. Polehampton   22-10-2013  - First version adapted from Andreas' script - SPCAL-83
#   E. Polehampton   31-10-2013  - update for new input file
#   E. Polehampton   01-11-2013  - Add beam area metadata
#   E. Polehampton   23-01-2014  - Update values (SPCAL-98)
#  Chris North   - 18/Feb/2014 - 1.0: First version (complete calculation)
#                                     Used for spire_cal_12_1
#                                 SPCAL-109
#===============================================================================
import os
scriptVersionString = "makeSCalPhotColorCorrBeam.py $Revision: 1.6 $"

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
#  e.g. <directory>/Phot/SCalPhotColorCorrK/SCalPhotColorCorrBeam_<version>.fits
# If outputCalDirTree is False, then calibration products are written to
#  dataDir e.g. <dataDir>/SCalPhotColorCorrBeam_<version>.fits
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
# Input parameters for colour correction
# range of alphas to compute colour corrections for
alphaK=[-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-0.5,0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
#alphaK=[-4,-1,0,2]
# range of beta and temp to calculate colour corrections for
betaK=[0.,0.5,1.,1.25,1.5,1.75,2.,2.5,3.]
tempK=range(3,300)

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
# Load SPIRE filter functions fron calibration directory tree

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
	#create Rsrf interpolation object
	interpRsrf = LinearInterpolator(spireFreq, rsrf.getRsrf(band))
	#make arrays for final object
	spireFiltOnly[band] = Double1d(nNu)
	#interpolate Rsrf to freq array
	spireFiltOnly[band][ixR] = interpRsrf(freq[ixR])

#-------------------------------------------------------------------------------
# Load SPIRE Beam profiles
beamProfsVersion = "3"
beamProfs = fitsReader("%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits"%(directory, beamProfsVersion))
spireEffFreq = {"PSW":beamProfs.meta['freqEffPsw'].double*1.e9,\
	"PMW":beamProfs.meta['freqEffPmw'].double*1.e9,\
	"PLW":beamProfs.meta['freqEffPlw'].double*1.e9}
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
#   * spireEffArea: Calculate the effective beam area for a given spectrum

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
# Calculate the effective beam area for a given spectrum
def spireEffArea(freq, transm, monoArea, BB=False, temp=20.0, beta=1.8, alpha=-1.0):
	"""
	========================================================================
	Calculate the effective beam area for a source of a given spectrum

	Inputs:
	  freq:       (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:     (array float) relative spectral response (RSRF) corresponding to freq
	                Note that this should *not* include the aperture efficiency
	  monoArea:   (array float) monochromatic beam solid angle corresponding
	                to frequencies in freq
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
	            (float) Beam area for given spectrum, in same units as monoArea

	Calculation:
	  Calculates the source spectrum (either modifies black body or power law)
          Multiplies the monochromatic beam area by RSRF and source spectrum
	  Integrates over frequency
	  Normalises by integral over frequency of RSRF and source spectrum

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2013/12/19  C. North  initial version

	"""	
	#
	# Calculate sky background model
	#
	if BB == 1:
		#print temp,beta,c,h,k
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha

	# Integrate monochromatic area over frequency, weighted by rsrf and fSky
	numInterp=LinearInterpolator(freq,transm * fSky * monoArea)
	denomInterp=LinearInterpolator(freq,transm * fSky)
	minFreq=min(freq)
	maxFreq=max(freq)
	integrator=TrapezoidalIntegrator(minFreq,maxFreq)
	numInteg=integrator.integrate(numInterp)
	denomInteg=integrator.integrate(denomInterp)
	effArea = numInteg / denomInteg

	return(effArea)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                           END OF FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#===============================================================================
#=====                  CALCULATE MONOCHROMATIC BEAM AREAS                 =====
#===============================================================================
#-------------------------------------------------------------------------------

#calculate monochromatic beam areas using full or simple beam treatment
beamMonoArea={}
print '\nCalculating monochromatic beam areas...'
for band in spireBands:
	#monochromatic beam areas
	beamMonoArea[band] = spireMonoAreas(freq, beamProfs, 
	  spireEffFreq[band], gamma, band)

# effBeamSr: * Set of tables per band for various source spectra
#            * Effective beam solid angle (in sr) for a source with a given spectrum
#            * These values are not stored in the calibration tree
#
# kBeam: * Set of tables per band for various source spectra
#        * Defined as beamAreaPipSr / effBeamSr
#        * Colour correction for beam solid angle
#        * These values are stores in the ColorCorrBeam tables
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# read pipeline values from beam profile calibration product

beamAreaPipSr  = {'PSW':beamProfs.meta['beamPipelinePswSr'],\
	'PMW':beamProfs.meta['beamPipelinePmwSr'],\
	'PLW':beamProfs.meta['beamPipelinePlwSr']}
beamAreaPipArc = {'PSW':beamProfs.meta['beamPipelinePswArc'],\
	'PMW':beamProfs.meta['beamPipelinePmwArc'],\
	'PLW':beamProfs.meta['beamPipelinePlwArc']}
#create dictionaries for filenames

#-----------------------------------------------------------------------
#=======================================================================
#=====                 CALCULATE PIPELINE PARAMETERS               =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Calculate pipeline colour correction parameters
print '\nGenerating new pipeline beam areas for ColorCorrBeam version %s'%version
for band in spireBands:
	#pipeline beam areas
	beamAreaPipSr[band]=spireEffArea(freq, spireFiltOnly[band], \
	  beamMonoArea[band], BB=False, alpha=-1)
	beamAreaPipArc[band]=beamAreaPipSr[band]/arcsec2Sr

#-----------------------------------------------------------------------
#=======================================================================
#=====               CALCULATE BEAM COLOR CORRECTIONS              =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#create product and set basic metadata
# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

kCorrBeam=herschel.spire.ia.dataset.PhotColorCorrBeam()
kCorrBeam.setDescription('SPIRE beam corrections with spectral index and temperature product')
kCorrBeam.meta["creator"].value   = scriptVersionString
kCorrBeam.meta["modelName"].value = "FM"
kCorrBeam.meta["creationDate"].value = FineTime(java.util.Date())
kCorrBeam.meta["startDate"].value = FineTime(startDate)
kCorrBeam.meta["endDate"].value   = FineTime(endDate)
kCorrBeam.meta["author"]  = herschel.ia.dataset.StringParameter(value="Chris North", description="Author of the data")
kCorrBeam.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="v3.0", description="Origin of the data")
kCorrBeam.setVersion(version)
kCorrBeam.setFormatVersion(formatVersion)

# make product contain raw effective beam areas
effBeamSr=herschel.spire.ia.dataset.PhotColorCorrBeam()
effBeamSr.setDescription('SPIRE beam solid angle with spectral index and temperature product')

#-----------------------------------------------------------------------
# Compute beam corrections for range of alpha and temp/beta
# Create tables for alpha arrays for KBeam
effBeamSr['alpha']=TableDataset(description='Beam Solid Angle (Spectral Index)')
kCorrBeam['alpha']=TableDataset(description='Beam Colour Correction (Spectral Index)')

#add alpha column
effBeamSr['alpha'].addColumn('alpha',Column(Float1d(alphaK)))
kCorrBeam['alpha'].addColumn('alpha',Column(Float1d(alphaK)))

print 'Calculating beam correction parameters over alpha...'
for band in spireBands:
	# add columns to tables
	effBeamSr['alpha'].addColumn(band,Column(Float1d(len(alphaK)),unit=SolidAngle.STERADIANS,description=''))
	kCorrBeam['alpha'].addColumn(band,Column(Float1d(len(alphaK))))

	for a in range(len(alphaK)):
		#effective area
		effBeamSr['alpha'][band].data[a]=spireEffArea(freq, spireFiltOnly[band],
		  beamMonoArea[band], BB=False, alpha=alphaK[a])
		#beam correction factor
		kCorrBeam['alpha'][band].data[a] = \
		  beamAreaPipSr[band] / effBeamSr['alpha'][band].data[a]

#-----------------------------------------------------------------------
print 'Calculating beam colour correction parameters over beta & temp...'
for b in range(len(betaK)):
	#create format version of beta text (beta=x.yz -> beta_x_yz)
	betaTxt='beta_%d_%d%d'%(int(betaK[b]),int(10*(betaK[b]%1)),int(10*((10*betaK[b])%1)))
	if verbose:
		print '  %s: beta=%f'%(betaTxt,betaK[b])
	#create tables
	effBeamSr[betaTxt]=TableDataset(description='Beam Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
	kCorrBeam[betaTxt]=TableDataset(description='Beam Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
	# add temp column
	effBeamSr[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
	kCorrBeam[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
	for band in spireBands:
		effBeamSr[betaTxt].addColumn(band,Column(Float1d(len(tempK)),unit=SolidAngle.STERADIANS,description=''))
		kCorrBeam[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
		for t in range(len(tempK)):
			#effective area
			effBeamSr[betaTxt][band].data[t]=spireEffArea(freq, spireFiltOnly[band],
			  beamMonoArea[band], BB=True, beta=betaK[b], temp=tempK[t])
			#beam correction factor
			kCorrBeam[betaTxt][band].data[t] = \
			  beamAreaPipSr[band] / effBeamSr[betaTxt][band].data[t]

#update ColorCorrBeam metadata
kCorrBeam.meta['beamPswSr']=DoubleParameter(beamAreaPipSr['PSW'],\
 unit=SolidAngle.STERADIANS,description='PSW beam area for spectral index alpha=-1 (as assumed in the pipeline')
kCorrBeam.meta['beamPmwSr']=DoubleParameter(beamAreaPipSr['PMW'],\
 unit=SolidAngle.STERADIANS,description='PMW beam area for spectral index alpha=-1 (as assumed in the pipeline')
kCorrBeam.meta['beamPlwSr']=DoubleParameter(beamAreaPipSr['PLW'],\
 unit=SolidAngle.STERADIANS,description='PLW beam area for spectral index alpha=-1 (as assumed in the pipeline')
kCorrBeam.meta['beamPswArc']=DoubleParameter(beamAreaPipArc['PSW'],\
 unit=SolidAngle.SQUARE_SECONDS_ARC,description='PSW beam area for spectral index alpha=-1 (as assumed in the pipeline')
kCorrBeam.meta['beamPmwArc']=DoubleParameter(beamAreaPipArc['PMW'],\
 unit=SolidAngle.SQUARE_SECONDS_ARC,description='PMW beam area for spectral index alpha=-1 (as assumed in the pipeline')
kCorrBeam.meta['beamPlwArc']=DoubleParameter(beamAreaPipArc['PLW'],\
 unit=SolidAngle.SQUARE_SECONDS_ARC,description='PLW beam area for spectral index alpha=-1 (as assumed in the pipeline')

# set FITS filename
if outputCalDirTree:
	kCorrBeamFits = java.io.File(r"%s//Phot//SCalPhotColorCorrBeam//SCalPhotColorCorrBeam_v%s.fits"%(directory, version))
else:
	kCorrBeamFits = java.io.File(r"%s//SCalPhotColorCorrBeam_v%s.fits"%(dataDir, version))
kCorrBeam.meta['fileName'] = herschel.ia.dataset.StringParameter(value=kCorrBeamFits.name,\
  description="Name of file when exported")

#-----------------------------------------------------------------------
# write to FITS file
fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(kCorrBeamFits.toString(), kCorrBeam)
