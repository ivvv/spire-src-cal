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
#    SCalPhotColorCorrK_point product
#    SCalPhotColorCorrK_extended product
#
#===============================================================================
# 
#  Herschel-SPIRE Colour Corrections
# 
#  This routine calculates the colour corrections for a range of source spectra,
#  and for both point sources and extended sources
# 
#  Input:
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles & metadata from SPIRE calibration tree
#    Beam colour corrections from SPIRE calibraiton tree
#    Neptune spectral index and measured area
#    Array of alphas and (beta,temperatures) from SPIRE calibration tree
#    Name and version of output file
# 
#  Output:
#    SCalPhotColorCorrK product
# 
#  Calculations:
#   1. Reads in beam profiles, beam corrections, RSRF and aperture efficiency
#      from calibration tree, including the values of gamma, effective beam area
#      and the pipeline beam areas
#   2. Calculates the point source colour corrections as a function of alpha:
#        Kpsrc = int{ RSRF(nu) * ApEff(nu) dnu} \
#                     int{ (nu/nu0)^alpha * RSRF(nu) * ApEff(nu) dnu}
#      And similarly for beta, temperatyre
#   3. Scales the monochromatic beam with frequency by exponent gamma and
#      calculates the monochromatic beam area
#        Beam_mono(r,nu) = Beam_full(r*(nu/nu_eff)^-gamma)
#           (meaning that the FWHM varies as (nu/nu_eff)^gamma
#        Omega_mono(nu) = int{Beam_mono(r,nu) * 2*pi*r dr}
#   4. Calculates the effective beam area from beam correction
#        Omega_eff(alpha) = Omega_eff(-1) / Kbeam(alpha)
#   5. Computes the extended source colour corrections
#        Kextd = Omega_eff(alpha) * int{ RSRF(nu) * ApEff(nu) dnu} \
#              int{ (nu/nu0)^alpha * Omega_mono (nu) * RSRF(nu) * ApEff(nu) dnu}
#      And similarly for greybody spectra for a range of betas and temperature
#
#  There are some functions defined (common between calibration products)
#   * spireMonoBeam: Calculate monochromatic beam profile & area at a given freq
#   * spireMonoAreas: Calculate monochromatic beam areas at range of frequencies
#   * hpXcalKcorr: Calculate K-correction parameters for given spectrum
#
#===============================================================================
# $Id: makeSCalPhotColorCorrK.py,v 1.9 2014/02/26 18:42:48 epoleham Exp $
# 
#  Edition History
#   E. Polehampton   22-10-2013  - First version adapted from Andreas' script - SPCAL-83
#   E. Polehampton   31-10-2013  - update for new input file
#   E. Polehampton   21-01-2014  - add table descriptions and update numbers (SPCAL-93)
#   Chris North      18-02-2014  - updated to use full calculation instead of reading csv files
#
#===============================================================================
import os
scriptVersionString = "makeSCalPhotColorCorrK.py $Revision: 1.9 $"

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
#  e.g. <directory>/Phot/SCalPhotColorCorrK/SCalPhotColorCorrK_point_<version>.fits
# If outputCalDirTree is False, then calibration products are written to
#  dataDir e.g. <dataDir>/SCalPhotColorCorrK_point_<version>.fits
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

#set verbose to print more information during processing
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
	# SPIRE aperture efficiency product from cal tree
	apertureEfficiencyVersion = "1"
	apertureEfficiency = fitsReader("%s//Phot//SCalPhotApertureEfficiency//SCalPhotApertureEfficiency_v%s.fits"%(directory, apertureEfficiencyVersion))
	if verbose:
		print 'Reading RSRF version %s from calibration directory tree'%(rsrfVersion)
		print 'Reading Aperture Efficiency version %s from calibration directory tree'%(apertureEfficiencyVersion)
else:
	#read RSRF and Aperture Efficiency from calibration tree
	rsrf=cal.getPhot().getProduct('Rsrf')
	rsrfVersion=rsrf.getVersion()
	apertureEfficiency = cal.getPhot().getProduct('ApertureEfficiency')
	apertureEfficiencyVersion=apertureEfficiency.getVersion()
	if verbose:
		print 'Reading RSRF version %s from calibration %s'%(rsrfVersion,cal.getVersion())
		print 'Reading Aperture Efficiency version %s from calibration %s'%(apertureEfficiencyVersion,cal.getVersion())

# Photometer RSRF
spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz
#indexes of freq in rsrf
ixR = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# Photometer Aperture Efficiency
spireApEffFreq = apertureEfficiency.getApertEffTable()["frequency"].data * 1e9 #comes in [GHz]
#indexes of freq in apEff
ixA = freq.where((freq>=MIN(spireApEffFreq)) & (freq<=MAX(spireApEffFreq)))

# spire RSRF only
spireFiltOnly={}
# spire RSRF * ApEff
spireFilt={}

#interpolate to freq array
for band in spireBands:
	#create Rsrf and ApEff interpolation objects
	interpRsrf = LinearInterpolator(spireFreq, rsrf.getRsrf(band))
	interpAp = LinearInterpolator(spireApEffFreq, apertureEfficiency.getApertEffTable()[band].data)
	#make arrays for final objects
	spireFiltOnly[band] = Double1d(nNu)
	spireFilt[band] = Double1d(nNu)
	#interpolate Rsrf to freq array
	spireFiltOnly[band][ixR] = interpRsrf(freq[ixR])
	#copy into Rsrf*ApEff array
	spireFilt[band] = spireFiltOnly[band].copy()
	spireFilt[band][ixA] = spireFilt[band][ixA] * interpAp(freq[ixA])

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
# Calculate K-correction parameters for given spectrum & source type
def hpXcalKcorr(freq0, freq, transm, BB=True, temp=20.0, beta=1.8, alpha=-1.0,
  ext=False, monoArea=None):
	"""
	================================================================================
	Calculation of the K-correction factor from isophotal flux to a monochromatic 
	flux-density at a given reference frequency (data to be multiplied!)
	This routine is needed by hpXcalColorCorr.py
	
	Inputs:
	  freq0:     (float) waveband reference frequency [Hz] for which monochromatic
	               flux-density is given
	  freq:      (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:    (array float) relative spectral response (RSRF) corresponding to freq
	  BB:        (boolean) spectral function to use for source spectrum:
	                'True': a modified black body
	                'False' a power-law with exponent alpha=-1
	                OPTIONAL. Default=False
	  temp:      (float) Dust/sky temperature [K] (only for modified black body)
			OPTIONAL. Deafult=20K; 
	                only for modified black body]
	  beta:      (float) Dust/sky spectral index (only for modified black body]
			OPTIONAL. Default=1.8
	  alpha:     (float) Exponent of power-law sky background model (only for
			power-law spectrum)
			OPTIONAL. Default=-1
	  ext:       (boolean) calculating for extended source
	                OPTIONAL. Default=False
	  monoArea:  (array float) Monochromatic Beam solid angle [Sr] corresponding
	                to freq.
	                OPTIONAL. Only required if ext=True

	Outputs:
	 (list)     [0]: K-correction factor
	            [1]: Sky emission at reference fequency (fSky0)
	
	Calculation:
	  Depending on the state of the input parameter BB, either the spectrum of a
	  Planck function multiplied by frequency to the power of beta, or a power-law
	  spectrum with spectral index alpha is calculated for all values in the vector 
	  frequ. In addition the same value is calculated at the discrete frequency 
	  freq0. Then the product of this value and the integral of the RSRF over all
	  frequencies, divided by the integral over all products of frequency and RSRF	
	  is calculated. Note that the integrals are coded as simple sums as the 
	  HIPE numeric integral doesn't allow too may discrete points and the RSRF
	  is sampled to quite some detail.

	  N.B. If ext=True, and monoBeam is in units sr, then this procedure outputs
	    K-correction factor in [Jy/sr per Jy/beam] and Sky emission in Jy/sr.
	    Units will change if a different input unit is used.
	
	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2012/04/18  L. Conversi  initial version in file hp_xcal_Total_v1.6.py.txt
	2012/08/22  B. Schulz    added comments, reformatted, renamed kCorr to hpXcalKcorr
	2012/08/24  B. Schulz    fixed defaults for temp and beta, changed inputs to frequency,
	                         updated header and comments, renamed inputFreq to freq
	                         and inputFilt to transm added powerlaw sky background
	2012/08/27  B. Schulz    brush -up on header and comments
	2013/03/28  B. Schulz    removed implicit limitation to fixed frequency interval
	2013/04/05  B. Schulz    implemented proper tabulated integration function
	2013/06/18  L. Conversi  implemented gamma parameter in the denominator
	                         (previsuly applie directly to RSRFs, i.e. in the numerator too)
	2013/12/03  C. North     corrected procedure include area where appropriate
			 	 NB: if ext=True , note output units
	
	================================================================================
	"""
	#
	# Calculate sky background model
	#
	# 1) As a modified Blackbody
	if BB == 1:
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
		fSky0 = 2*h * freq0**3 / c**2 / (EXP(h*freq0/k/temp) - 1.) * freq0**beta
	#
	# 2) As a Power-Law
	else:
		fSky  = freq**alpha
		fSky0 = freq0**alpha
	#
	# Using the K-correction as defined in th OM (eq. 5.5, 5.15 & 5.16)

	if ext == True:
		area = monoArea
	else:
		#don't use area for point sources
		area = Float1d(len(freq))
		area[:] = 1.0

        # monoArea is the monochromatic beam solid angle at effFreq

	# integrate over frequency
	numInterp=LinearInterpolator(freq,transm)
	denomInterp=LinearInterpolator(freq,transm * fSky * area)
	minFreq=min(freq)
	maxFreq=max(freq)

	integrator=TrapezoidalIntegrator(minFreq,maxFreq)
	numInteg = integrator.integrate(numInterp)
	denomInteg = integrator.integrate(denomInterp)

	kWave = fSky0 * numInteg / denomInteg

	#
	# Return the result as a 2-element array of K-correction and flux at freq0
	return (Double1d([kWave, fSky0]))

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
beamAreaPipSr  = {}
beamAreaPipArc  = {}
print '\nCalculating monochromatic beam areas...'
for band in spireBands:
	#monochromatic beam areas
	beamMonoArea[band] = spireMonoAreas(freq, beamProfs, 
	  spireEffFreq[band], gamma, band)

#-------------------------------------------------------------------------------
#===============================================================================
#=====                      CALCULATE COLOR CORRECTIONS                    =====
#===============================================================================
#-------------------------------------------------------------------------------
# Description of parameters:
# K4P: * Scalar value per band
#      * Converts from broadband RSRF-weighted flux density (Jy/beam) to
#        monochromatic flux density (Jy/beam) at ref frequency a for point source
#        with alpha=-1 power law spectrum (nu*Fnu = const)
#      * This is applied in pipeline for the point source (PxWpsrc) L2 products
#      * These are stored in FluxConv metadata
#
# beamAreaPipSr: * Scalar value per band.
#                * Effective beam area (in sr) for a source with alpha=-1 spectrum
#                * These are stored in RadialCorrBeam and ColorCorrBeam metadata
#
# beamAreaPipArc: * Scalar value per band.
#                 * Effective beam area (in arcsec^2) for a source with alpha=-1 spectrum
#                 * These are stored in RadialCorrBeam and ColorCorrBeam metadata
#
# K4E_Tot: * Scalar value per band.
#          * Converts from broadband RSRF-weighted flux density (Jy/beam) to
#            monochromatic surface brightness (Jy/sr) at reference frequency for an
#            extended source with alpha=-1 spectrum (i.e. nu*Fnu = const)
#          * This is essentially what is applied in pipeline for the extended
#            source (PxWextd) L2 products. In actualitly, point source flux 
#            (with K4P applied) is multiplied by K4E * beamAreaPipSr / K4P.
#          * These values are not stored in the calibration tree
#
# K4E: * Scalar value per band.
#      * Defined as K4E_Tot * beamAreaPipSr.
#      * Converts from broadband RSRF-weighted flux density (Jy/beam) to
#        monochromatic flux density (Jy/beam) at reference frequency for an
#        extended source with alpha=-1 spectrum (i.e. nu*Fnu = const)
#      * These are stored in FluxConv metadata
#
# KPtoE: * Scalar value per band.
#        * Defined as K4E_Tot / K4P.
#        * Converts from monochromatic flux density (Jy/beam) of a point source
#          with alpha=-1 spectrum to the monochromatic surface brightness (Jy/sr)
#          of an extended source with alpha=-1 spectrum.
#        * These values are not stored in the calibration tree
#
# KcP: * Set of tables per band for various source spectra
#      * Colour correction to convert from monochromatic flux density (Jy/beam)
#        of a point source with alpha=-1 spectrum to the monochromatic flux density
#        (Jy/beam) of a point source with a given spectrum.
#      * These values are stored in the ColorCorrK_point tables
#
# KcE: * Set of tables per band for various source spectra
#      * Colour correction to convert from monochromatic flux density (Jy/beam)
#        of an extended source with alpha=-1 spectrum to the monochromatic flux 
#        density (Jy/beam) of an extended source with a given spectrum.
#      * These values are stored in the ColorCorrK_extended tables
#
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
# Compute colour corrections for range of alpha and temp/beta
#set up dictionary for pipeline values
k4P = {}
k4E_Tot = {}
k4E = {}
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
	#pipeline point source K4P conversion
	k4P[band] = hpXcalKcorr(spireRefFreq[band], freq, spireFilt[band], \
	  BB=False, ext=False)[0]
	#pipeline extended source K4E conversion
	k4E_Tot[band] = hpXcalKcorr(spireRefFreq[band], freq, spireFilt[band],
	  BB=False, ext=True, monoArea=beamMonoArea[band])[0]
	k4E[band] = k4E_Tot[band] * beamAreaPipSr[band]

#-----------------------------------------------------------------------
#=======================================================================
#=====           CALCULATE POINT SOURCE COLOR CORRECTIONS          =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#create product and set basic metadata
# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

kCorrPsrc = herschel.spire.ia.dataset.PhotColorCorrK()
kCorrPsrc.description = "SPIRE Color Corrections product for point sources"
kCorrPsrc.meta["creator"].value   = scriptVersionString
kCorrPsrc.meta["modelName"].value = "FM"
kCorrPsrc.meta["creationDate"].value = FineTime(java.util.Date())
kCorrPsrc.meta["startDate"].value = FineTime(startDate)
kCorrPsrc.meta["endDate"].value   = FineTime(endDate)
kCorrPsrc.meta["author"]  = herschel.ia.dataset.StringParameter(value="Chris North", description="Author of the data")
kCorrPsrc.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="v3.0", description="Origin of the data")
kCorrPsrc.meta["dependency"].value = "correctionSourceType"
kCorrPsrc.meta["correctionSourceType"] = StringParameter(value="point", description="")
kCorrPsrc.setVersion(version)
kCorrPsrc.setFormatVersion(formatVersion)

#make dummy product for point source conversion
#kConvPsrc=herschel.spire.ia.dataset.PhotColorCorrK()
#kConvPsrc.setDescription('Spire color conversions with spectral index and temperature product')

print '\nCalculating point source colour correction parameters over alpha...'
# Create tables for alpha arrays for K (point source)
kCorrPsrc.setAlphaValues(Double1d(alphaK))

# Add description to table
kCorrPsrc['alpha'].description = "Point Source Colour Correction (Spectral Index)"

#kConvPsrc['alpha']=TableDataset(description='Point Source Conversion (Spectral Index)')
# add alpha column
#kConvPsrc['alpha'].addColumn('alpha',Column(Float1d(alphaK)))

for band in spireBands:
	# add band columns to tables
	alphaData = Double1d(len(alphaK))
	#kConvPsrc['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
	for a in range(len(alphaK)):
		#point source conversion for current alpha
		kConvPsrc=hpXcalKcorr(spireRefFreq[band],\
		   freq, spireFilt[band], BB=False, alpha=alphaK[a])[0]
		#point source colour correction for current alpha
		alphaData[a] = kConvPsrc/k4P[band]
	kCorrPsrc.setCorrVals(band, alphaData)


#-----------------------------------------------------------------------
print 'Calculating point source colour correction parameters over beta & temp...'
for b in range(len(betaK)):
	#create format version of beta text (beta=x.yz -> beta_x_yz)
	betaTxt='beta_%d_%d%d'%(int(betaK[b]),int(10*(betaK[b]%1)),int(10*((10*betaK[b])%1)))
	if verbose:
		print '  %s: beta=%f'%(betaTxt,betaK[b])
	# Create tables for beta arrays for K (point source)
	#kConvPsrc[betaTxt]=TableDataset(description='Point Source Conversion (Modified Black Body, beta=%.2f)'%betaK[b])
	# add temp column
	kCorrPsrc.setTempBetaVals(betaK[b], Double1d(tempK))
	#kConvPsrc[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
	for band in spireBands:
		betaData = Double1d(len(tempK))
		#kConvPsrc[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
		for t in range(len(tempK)):
			#point source conversion for current beta,temp
			kConvPsrc=hpXcalKcorr(spireRefFreq[band],\
			  freq, spireFilt[band], BB=True, beta=betaK[b], temp=tempK[t])[0]
			#point source colour correction for current beta,temp
			betaData[t] = kConvPsrc/k4P[band]
		kCorrPsrc.setTempBetaCorrVals(betaK[b], band, betaData)
	kCorrPsrc.getTempBetaCorrTable(betaK[b]).description = "Point Source Colour Correction (Modified Black Body, beta=%.2f)"%(betaK[b])

#-----------------------------------------------------------------------
#=======================================================================
#=====         CALCULATE EXTENDED SOURCE COLOR CORRECTIONS         =====
#=======================================================================
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

# update basic metadata
kCorrExtd = herschel.spire.ia.dataset.PhotColorCorrK() 
kCorrExtd.description = "SPIRE Color Corrections product for extended sources"
kCorrExtd.meta["creator"].value   = scriptVersionString
kCorrExtd.meta["modelName"].value = "FM"
kCorrExtd.meta["creationDate"].value = FineTime(java.util.Date())
kCorrExtd.meta["startDate"].value = FineTime(startDate)
kCorrExtd.meta["endDate"].value   = FineTime(endDate)
kCorrExtd.meta["author"]  = herschel.ia.dataset.StringParameter(value="Chris North", description="Author of the data")
kCorrExtd.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="v3.0", description="Origin of the data")
kCorrExtd.meta["dependency"].value = "correctionSourceType"
kCorrExtd.meta["correctionSourceType"] = StringParameter(value="extended", description="")
kCorrExtd.setVersion(version)
kCorrExtd.setFormatVersion(formatVersion)

#create dummy for extended conversion parameter
#kConvExtd = herschel.spire.ia.dataset.PhotColorCorrK()
#kConvExtd.setDescription('Spire conversions with spectral index and temperature product')
#-----------------------------------------------------------------------
# Compute extended source colour corrections for range of alpha and temp/beta
# Create tables for alpha arrays for K (extended source)
kCorrExtd.setAlphaValues(Double1d(alphaK))

# Add description to table
kCorrExtd['alpha'].description = "Extended Source Colour Correction (Spectral Index)"

#kConvExtd['alpha']=TableDataset(description='Extended Source Conversion (Spectral Index)')
# add alpha column
#kConvExtd['alpha'].addColumn('alpha',Column(Float1d(alphaK)))

#-----------------------------------------------------------------------
print '\nCalculating extended source colour correction parameters over alpha...'
for band in spireBands:
	# add band columns to tables
	alphaData = Double1d(len(alphaK))
	#kConvExtd['alpha'].addColumn(band,Column(Float1d(len(alphaK))))
	for a in range(len(alphaK)):
		#total extended source conversion for current alpha
		k4EaTot_x=hpXcalKcorr(spireRefFreq[band], freq,\
		 spireFilt[band], BB=False, alpha=alphaK[a],\
		 ext=True, monoArea=beamMonoArea[band])[0]
		#effective beam area
		effBeamSr = beamAreaPipSr[band] / kBeam['alpha'][band].data[a]
		#total conversion
		#kConvExtd['alpha'][band].data[a] = \
		#  k4EaTot_x * effBeamSr
		#extended source colour correction for current alpha
		alphaData[a] = k4EaTot_x / k4E_Tot[band]
	kCorrExtd.setCorrVals(band, alphaData)

#-----------------------------------------------------------------------
print 'Calculating extended source colour correction parameters over beta & temp...'
for b in range(len(betaK)):
	#create format version of beta text (beta=x.yz -> beta_x_yz)
	betaTxt='beta_%d_%d%d'%(int(betaK[b]),int(10*(betaK[b]%1)),int(10*((10*betaK[b])%1)))
	if verbose:
		print '  %s: beta=%f'%(betaTxt,betaK[b])
	# Create tables for beta for K (extended source)
	#kCorrExtd[betaTxt]=TableDataset(description='Extended Source Colour Correction (Modified Black Body, beta=%.2f)'%betaK[b])
	#kConvExtd[betaTxt]=TableDataset(description='Extended Source Conversion (Modified Black Body, beta=%.2f)'%betaK[b])
	# add temp column
	kCorrExtd.setTempBetaVals(betaK[b], Double1d(tempK))
	#kCorrExtd[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
	#kConvExtd[betaTxt].addColumn('Temperature',Column(Float1d(tempK),unit=Temperature.KELVIN))
	for band in spireBands:
		betaData = Double1d(len(tempK))
		#kConvExtd[betaTxt].addColumn(band,Column(Float1d(len(tempK))))
		for t in range(len(tempK)):
			#extended source conversion for current beta,temp
			k4EbTot_x=hpXcalKcorr(spireRefFreq[band],\
			  freq, spireFilt[band], BB=True, beta=betaK[b], temp=tempK[t],\
			  ext=True, monoArea=beamMonoArea[band])[0]
			#effective beam area
			effBeamSr = beamAreaPipSr[band] / kBeam[betaTxt][band].data[t]
			#kConvExtd[betaTxt][band].data[t] = \
			#  k4EbTot_x * effBeamSr
			#extended source colour correction for current beta,temp
			betaData[t] = k4EbTot_x / k4E_Tot[band]
		kCorrExtd.setTempBetaCorrVals(betaK[b], band, betaData)
	kCorrExtd.getTempBetaCorrTable(betaK[b]).description = "Extended Source Colour Correction (Modified Black Body, beta=%.2f)"%(betaK[b])		



# set FITS filenames
if outputCalDirTree:
	kCorrPsrcFits = java.io.File(r"%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_point_v%s.fits"%(directory, version))
	kCorrExtdFits = java.io.File(r"%s//Phot//SCalPhotColorCorrK//SCalPhotColorCorrK_extended_v%s.fits"%(directory, version))
else:
	kCorrPsrcFits = java.io.File(r"%s//SCalPhotColorCorrK_point_v%s.fits"%(dataDir, version))
	kCorrExtdFits = java.io.File(r"%s//SCalPhotColorCorrK_extended_v%s.fits"%(dataDir, version))
kCorrPsrc.meta['fileName'] = herschel.ia.dataset.StringParameter(value=kCorrPsrcFits.name,\
  description="Name of file when exported")
kCorrExtd.meta['fileName'] = herschel.ia.dataset.StringParameter(value=kCorrExtdFits.name,\
  description="Name of file when exported")

#-----------------------------------------------------------------------
# write to FITS files
fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(kCorrPsrcFits.toString(), kCorrPsrc)

fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(kCorrExtdFits.toString(), kCorrExtd)
