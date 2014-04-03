#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
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
# $Id: makeSCalPhotColorCorrHfi.py,v 1.3 2014/02/26 18:06:39 epoleham Exp $
#
#===============================================================================
# 
#  Planck-HFI / Herschel-SPIRE Colour Correction factors 
# 
#  This routine calculates a table of colour correction factors for a range of
#  realistic celestial background spectra. The factors transfer standard 
#  monochromatic fluxdensities from Planck-HFI all-sky maps at 857 and 545 GHz
#  which are quoted for a nu*F_nu = const. spectrum into Herschel-SPIRE 
#  monochromatic flux densities with the same definition but reference 
#  wavelengths at 250, 350, and 500 micron. It is assumed that the actual 
#  celestial background is described by the Planck-Law multiplied by frequency 
#  to the power of beta which is set to 1.8. In addition the flux ratios as 
#  measured by HFI between the 545 and 857 Gz filters is given such that the 
#  colour corrections can be determined directly as a function of the flux 
#  ratios in those maps.
# 
#  Input:
#    HFI RIMO files (containing HFI RSRF profiles), downloadable at:
#        http://pla.esac.esa.int/pla/aio/product-action?DOCUMENT.DOCUMENT_ID=HFI_RIMO_R1.10.fits
#    RSRF and Aperture efficiency profiles from SPIRE calibration tree
#    Beam profiles and beam model parameters from SPIRE calibration tree
#     - These can be from a calibration tree structure or a SpireCal object
#     - Set inputCalDirTree and outputCalDirTree accordingly
#    Flux conversion parameters from SPIRE calibration tree
#    Exponent for modified Blackbody to model dust emission background colour
#    Array of temperatures from SPIRE calibration tree
#    Name and version of output file
# 
#  Output:
#    SCalPhotColorCorrHfi product
# 
#  Calculations:
#   K-Correction factors to go from HFI to SPIRE maps. The colour correction 
#   computed here:
#   1. starts from HFI maps computed under the assumption of having extended 
#      source with alpha=-1 spectrum
#   2. "converts" them to sky flux density assuming a grey body spectrum
#   3. "converts" them again to SPIRE wavebands using the same grey body 
#      assumption
#   4. finally changes them to the SPIRE extended emission calibration 
#      assuming an alpha=-1 spectrum
#  
#   In other terms: 
#   kHFItoSPIRE = 
#      1/K4E_HFI * K4E_MOD_HFI / K4E_MOD_SPIRE * K4E_SPIRE * fSky_SPIRE / fSky_HFI
#  
#   We also calculate the ratio of the Planck maps as measured by the pipeline,
#   i.e. under the assumption alpha=-1, but assuming a grey body of temperature T 
#   and spectral index beta. Hence:
#   R_HFI = fPipe_545 / fPipe_857 =
#         = (fSky_GB_545 / k545 * k4E_545) / (fSky_GB_857 / k857 * k4E_857) 
#  
#   To do this we require the monochromatic beam area as a function of frequency.
#   The SPIRE beam profiles area published based on observations of Neptune. An
#   effective frequency is defined as the frequency at which the monochromatic
#   beam area matches the broadband RSRF-weighted beam area measured on Neptune.
#   The beam profile is scaled as nu*gamma, and the resulting area computed
# 
# 
#  Edition History
#  Luca Conversi   - 13/Jul/2012 - 0.1: First test version
#  Bernhard Schulz - 18/Aug/2012 - 0.2: Separated out and added plots
#  Bernhard Schulz - 22/Aug/2012 - 0.3: Added comments, reformatted, file output added,
#                                       renamed kCorr to hpXcalKcorr
#  Bernhard Schulz - 24/Aug/2012 - 0.4: Addded point source RSRF to plot and calculate 
#                                       K4 factors for point sources
#  Bernhard Schulz - 27/Aug/2012 - 0.5: Cleaned up code and added more comments, added 
#                                       point source and extended source colour correction 
#                                       factors for standard spectrum as metadata
#  Bernhard Schulz - 28/Aug/2012 - 1.0: Renamed ratioPMWoverPLW to kPMWtoPLW and included 
#                                       into output file as column PMWtoPLW
#  Luca Conversi   - 16/Oct/2012 - 1.1: Change column names and order;
#                                       computing k857toPSW instead of kPMWtoPSW
#  Luca Conversi   - 21/Mar/2013 - 2.0: Adding the possibility of constructing a
#                                       colour correction table with T fixed and variable beta
#                                       Other minor changes (e.g. toggles plots on/off)
#  Bernhard Schulz - 05/Apr/2013 - 2.1: Change to tabular integration, 
#                                       added gamma parameter to header.
#  Bernhard Schulz - 12/Apr/2013 - 2.2: Common frequency interpolation and aperture 
#                                       efficiencies included
#  Luca Conversi   - 30/May/2013 - 2.3: Making file indipendent from external librabry
#                                       Changing script to read HFI RSRF from "RIMO" file
#                                       (Table of the performance and instrumental characteristics
#                                       as downloaded from the Planck legacy archive)
#                                       Cleaning of filenames, paths, etc.
#  Luca Conversi   - 18/Jun/2013 - 2.4: Update Ks computation following C. North inputs:
#                                       now gamma is kept separatated from SPIRE RSRFs definition
#                                       and included in the denominator of the Ks calculations
#                                       Correct typo in table metadata when fixBeta = False
#  Ed Polehampton  - 23/Oct/2013 - 2.5  Adapt to produce the new SCalPhotColorCorrHfi calibration
#                                       product (SPCAL-77)
#                                       Removed K4 parameters from metadata
#  Chris North     - 15/Jan/2014 - 2.6: Implemented full beam model treatment
#                                       Reads in RadialCorrBeam calibration product
#                                       Reads in K4P and K4E from FluxConv calibration product
#                                       Added functions to compute effective beam areas
#                                       Changed frequency range to include all of PSW band
#                                       Changed temperature range to match ColorCorrK products
#                                       Updated K computation
#                                       Changed integration method to TrapezoidalIntegrator
#                                       Changed RSRF/Aperture Efficiency interpolation to Linear
#  Chris North     - 21/Feb/2014 - 2.7: Added option to read inputs from calProduct
#                                       BUG FIX: Corrected typo in 857-PMW conversion
#  Chris North     - 25/Feb/2014 - 2.8: BUG FIX: Uses KmonE not K4E in calculations (SPCAL-109)
#                                       
#===============================================================================

scriptVersionString = "makeSCalPhotColorCorrHfi.py $Revision: 1.3 $"

directory = "..//..//..//..//..//..//data//spire//cal//SCal"
dataDir = "//disks//winchester2//calibration_data//"
#LOCAL VERSION
#directory = Configuration.getProperty('var.hcss.workdir')
#dataDir = Configuration.getProperty('var.hcss.workdir')

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
# Input parameters

# Colour correction table version
version = "3"

# set format version and date format
formatVersion = "1.0"
df  = java.text.SimpleDateFormat("yyyy.MM.dd/HH:mm:ss/z")


# Choose if you want to keep constant either:
# - grey body spectral index beta: fixBeta = True [default]
# - or grey body temperature: fixBeta = False
fixBeta = True

# Grey body spectral index used for simulated background dust emission spectrum
# Only used if fixBeta = True
beta0 = 1.8  

# Grey body temeprature (in K) used for simulated background dust emission spectrum
# Only used if fixBeta = False
temp0 = 20.0

# HFI RIMO file
hfiRIMOFileName = 'HFI_RIMO_R1.10.fits'
hfiRIMOFile = '%s//%s'%(dataDir,hfiRIMOFileName)

if inputCalDirTree:
	# SPIRE Photometer RSRF calibration product from cal directory tree
	rsrfVersion = "2"
	rsrf = fitsReader("%s//Phot//SCalPhotRsrf//SCalPhotRsrf_v%s.fits"%(directory, rsrfVersion))
	# SPIRE aperture efficiency product from cal directory tree
	apertureEfficiencyVersion = "1"
	apertureEfficiency = fitsReader("%s//Phot//SCalPhotApertureEfficiency//SCalPhotApertureEfficiency_v%s.fits"%(directory, apertureEfficiencyVersion))
 
	# SPIRE Photometer radial beam profiles from cal directory tree
	beamProfsVersion = "3"
	beamProfs = fitsReader("%s//Phot//SCalPhotRadialCorrBeam//SCalPhotRadialCorrBeam_v%s.fits"%(directory, beamProfsVersion))

	# SPIRE Photometer FluxConv from cal directory tree (for K4P and K4E parameters)
	#>>>CHECKL is this the right way to read in one FluxConv product<<<
	fluxConvVersion = "9"
	fluxConv = fitsReader("%s//Phot//SCalPhotFluxConv//SCalPhotFluxConv_nominal_20090819_v%s.fits"%(directory, fluxConvVersion))
else:
	rsrf=cal.getPhot().getProduct('Rsrf')
	rsrfVersion=rsrf.getVersion()
	apertureEfficiency = cal.getPhot().getProduct('ApertureEfficiency')
	apertureEfficiencyVersion=apertureEfficiency.getVersion()
	beamProfs = cal.getPhot().getProduct('RadialCorrBeam')
	beamProfsVersion=beamProfs.getVersion()
	fluxConv=cal.getPhot().getProduct('FluxConvList').refs[0]
	fluxConvVersion=fluxConv.getVersion()
	if verbose:
		print 'Reading RSRF version %s from calibration %s'%(rsrfVersion,cal.getVersion())
		print 'Reading Aperture Efficiency version %s from calibration %s'%(apertureEfficiencyVersion,cal.getVersion())
		print 'Reading Beam Radial Profiles version %s from calibration %s'%(beamProfsVersion,cal.getVersion())
		print 'Reading Flux Conversion version %s from calibration %s'%(fluxConvVersion,cal.getVersion())

# Choose if you want to plot the computed RSRFs and colour correction parameters
plot = True

#-------------------------------------------------------------------------------
# Loading physical constants
from herschel.share.unit import *
h = Constant.H_PLANCK.value
k = Constant.K_BOLTZMANN.value
c = Constant.SPEED_OF_LIGHT.value

# Frequency raster of common frequency grid
deltaNu = 0.1e9		# 0.1 GHz
nuMin   = 150.e9
nuMax   = 1800.e9
nNu     = FIX((nuMax-nuMin)/deltaNu)
freq    = Double1d(range(nNu)) * deltaNu + nuMin

# Define the temperatures range and number of rows in output table
# read temperature range from ColorCorrK
if fixBeta == True:
	Tmin = 5.		# Min. modified blackbody temperature
	Tmax = 300.		# Max. modified blackbody temperature
	ndiv = 500		# Number of elements in temperature vector

# Define the spectral index range and number of rows in output table
if fixBeta == False:
	Bmin = 1.		# Min. modified blackbody temperature
	Bmax = 3.		# Max. modified blackbody temperature
	ndiv = 200		# Number of elements in temperature vector

# Beam model parameters read from RadialCorrBeam
# Exponent of powerlaw describing FWHM dependence on frequency
# FWHM ~ frequ**gamma
gamma = beamProfs.meta['gamma'].double
#pipeline beam areas (in sr)
arcsec2Sr = (Math.PI/(60.*60.*180))**2
beamAreaPipPSW = beamProfs.meta['beamPipelinePswSr'].double
beamAreaPipPMW = beamProfs.meta['beamPipelinePmwSr'].double
beamAreaPipPLW = beamProfs.meta['beamPipelinePlwSr'].double
# effective frequencies
# monochromatic beam area = Neptune beam area at nu=effective frequency
spireEffFreqPSW = beamProfs.meta['freqEffPsw'].double*1.e9
spireEffFreqPMW = beamProfs.meta['freqEffPmw'].double*1.e9
spireEffFreqPLW = beamProfs.meta['freqEffPlw'].double*1.e9


# Three SPIRE filter reference frequencies for PSW, PMW, PLW respectively
spireRefFreq = c/Double1d([250.,350.,500.])*1e6 

# Two HFI filter reference frequencies for the 857 and 545 GHz filters respectively
hfiRefFreq = Double1d([857.,545.])*1e9       


#-------------------------------------------------------------------------------
# Load SPIRE filter functions

# Photometer RSRF
spireFreq   = rsrf['rsrf']['frequency'].data*1e9  # Frequency in Hz

ix = freq.where((freq>=MIN(spireFreq)) & (freq<=MAX(spireFreq)))

# Get RSRF for normal point sources
# Interpolate to common frequency grid
interpPLW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['plw'].data)
interpPMW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['pmw'].data)
interpPSW = CubicSplineInterpolator(spireFreq, rsrf['rsrf']['psw'].data)

spireFiltOnlyPLW = Double1d(nNu)
spireFiltOnlyPMW = Double1d(nNu)
spireFiltOnlyPSW = Double1d(nNu)

spireFiltOnlyPLW[ix] = interpPLW(freq[ix])
spireFiltOnlyPMW[ix] = interpPMW(freq[ix])
spireFiltOnlyPSW[ix] = interpPSW(freq[ix])

# Aperture efficiency table

spireApEffFreq = apertureEfficiency['frequency']['frequency'].data * 1e9 #comes in [GHz]
spireApEffPsw  = apertureEfficiency['frequency']["PSW"].data
spireApEffPmw  = apertureEfficiency['frequency']["PMW"].data
spireApEffPlw  = apertureEfficiency['frequency']["PLW"].data

# Fold in and interpolate aperture efficiency
ix = freq.where((freq>=MIN(spireApEffFreq)) & (freq<=MAX(spireApEffFreq)))
interpPLW = CubicSplineInterpolator(spireApEffFreq, spireApEffPlw)
interpPMW = CubicSplineInterpolator(spireApEffFreq, spireApEffPmw)
interpPSW = CubicSplineInterpolator(spireApEffFreq, spireApEffPsw)

spireFiltPLW = Double1d(nNu)
spireFiltPMW = Double1d(nNu)
spireFiltPSW = Double1d(nNu)

spireFiltPLW[ix] = interpPLW(freq[ix]) * spireFiltOnlyPLW[ix]
spireFiltPMW[ix] = interpPMW(freq[ix]) * spireFiltOnlyPMW[ix]
spireFiltPSW[ix] = interpPSW(freq[ix]) * spireFiltOnlyPSW[ix]


#-------------------------------------------------------------------------------
# Load HFI filter functions
hfiRIMO = fitsReader(file = hfiRIMOFile)

# Read and parse HFI filter files into tables 
hfiFreq545 = 100. * c * hfiRIMO['BANDPASS_F545']['WAVENUMBER'].data[1:-1]
hfiFreq857 = 100. * c * hfiRIMO['BANDPASS_F857']['WAVENUMBER'].data[1:-1]

hfiTrans545 = Double1d(hfiRIMO['BANDPASS_F545']['TRANSMISSION'].data[1:-1])
hfiTrans857 = Double1d(hfiRIMO['BANDPASS_F857']['TRANSMISSION'].data[1:-1])


# Excluding data points that are not monotonically increasing for HFI-545
diff = Double1d(len(hfiFreq545)-1)
for i in range(len(hfiFreq545)-1):
	diff[i] = hfiFreq545[i+1]-hfiFreq545[i]

ind = diff.where(diff <= 0.).toInt1d()
for i in ind:
	hfiFreq545.delete(i+1,1)
	hfiTrans545.delete(i+1,1)


# Excluding data points that are not monotonically increasing for HFI-545
diff = Double1d(len(hfiFreq857)-1)
for i in range(len(hfiFreq857)-1):
	diff[i] = hfiFreq857[i+1] - hfiFreq857[i]

ind = diff.where(diff <= 0.).toInt1d()
for i in ind:
	hfiFreq857.delete(i+1,1)
	hfiTrans857.delete(i+1,1)


# Interpolate HFI transmissions to common frequency grid 
hfiFilt545 = Double1d(nNu)
hfiFilt857 = Double1d(nNu)

ix545     = freq.where((freq>=MIN(hfiFreq545)) & (freq<=MAX(hfiFreq545)))
interp545 = CubicSplineInterpolator(hfiFreq545, hfiTrans545)
hfiFilt545[ix545] = interp545(freq[ix545])

ix857     = freq.where((freq>=MIN(hfiFreq857)) & (freq<=MAX(hfiFreq857)))
interp857 = CubicSplineInterpolator(hfiFreq857, hfiTrans857)
hfiFilt857[ix857] = interp857(freq[ix857])

# Normalize transmissions
hfiFilt545 = hfiFilt545 / MAX(hfiFilt545)
hfiFilt857 = hfiFilt857 / MAX(hfiFilt857)

#-------------------------------------------------------------------------------
# Plot relative spectral response functions

if plot:
	# Make linear plot of all filter bands including those for point sources
	p = PlotXY()
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt857, name='hfiFilt857'))
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt545, name='hfiFilt545'))
	#
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPSW * (freq/spireRefFreq[0])**(2*gamma), name='spireXFiltPSW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPMW * (freq/spireRefFreq[1])**(2*gamma), name='spireXFiltPMW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPLW * (freq/spireRefFreq[2])**(2*gamma), name='spireXFiltPLW'))
	#
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPSW, name='spireFiltPSW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPMW, name='spireFiltPMW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPLW, name='spireFiltPLW'))
	#
	p.xaxis.range = [100,800]
	p.yaxis.range = [-0.1,1.1]
	p.xaxis.titleText = "Wavelength [micron]"
	p.yaxis.titleText = "Relative Spectral Response"
	p.legend.visible = 1
	#
	# Make log plot in y direction of all extended source filter bands
	p = PlotXY()
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt857, name='hfiFilt857'))
	p.addLayer(LayerXY(c/freq*1e6, hfiFilt545, name='hfiFilt545'))
	#
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPSW * (freq/spireRefFreq[0])**(2*gamma), name='spireXFiltPSW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPMW * (freq/spireRefFreq[1])**(2*gamma), name='spireXFiltPMW'))
	p.addLayer(LayerXY(c/freq*1e6, spireFiltPLW * (freq/spireRefFreq[2])**(2*gamma), name='spireXFiltPLW'))
	#
	p.xaxis.range = [0,2000]
	p.yaxis.range = [1e-8,1.1]
	p.xaxis.titleText = "Wavelength [micron]"
	p.yaxis.titleText = "Relative Spectral Response"
	p.yaxis.type = Axis.LOG
	p.legend.visible = 1

#-------------------------------------------------------------------------------
# spireMonoBeam function definition
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
# spireMonoAreas function definition
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
# spireEffArea function definition
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
# hpXcalKcorr function definition
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

	  If ext=False then this procedure outputs dimensionless conversion factor
	  If ext=True then this procedure outputs K-correction factor to convert
	    to surface brightness (i.e. if monoBeam is in sr, conversion will be
	    [Jy/sr per Jy/beam]). The units of the Sky emission are unchanged.
	
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
                                 removed gamma keyword
                                 added ext keyword to specifiy if source is extended
                                 added monoArea keyword, used if ext=True
                                 output is now a Double1d, not list
                                 NB: if ext=False, still outputs conversion to flux density
			 	 NB: if ext=True , outputs conversion to surface brightness
                                     this needs to be divided by effective beam area for
                                     spectrum in question
	
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

	# integrate integrands of numerator and denominator over frequency
	numInterp=LinearInterpolator(freq,transm)
	denomInterp=LinearInterpolator(freq,transm * fSky * area)

	# perform integration
	minFreq=min(freq)
	maxFreq=max(freq)
	integrator=TrapezoidalIntegrator(minFreq,maxFreq)

	numInteg = integrator.integrate(numInterp)
	denomInteg = integrator.integrate(denomInterp)

	# calculate conversion factor
	kWave = fSky0 * numInteg / denomInteg

	# Return the result as a 2-element array of K-correction and flux at freq0
	return (Double1d([kWave, fSky0]))

#-------------------------------------------------------------------------------
print 'Calculating PSW monochromatic beam areas...'
beamMonoAreaPSW = spireMonoAreas(freq, beamProfs, spireEffFreqPSW, gamma, 'PSW')
print 'Calculating PMW monochromatic beam areas...'
beamMonoAreaPMW = spireMonoAreas(freq, beamProfs, spireEffFreqPMW, gamma, 'PMW')
print 'Calculating PLW monochromatic beam areas...'
beamMonoAreaPLW = spireMonoAreas(freq, beamProfs, spireEffFreqPLW, gamma, 'PLW')

#-------------------------------------------------------------------------------
# Compute K-correction factors for SPIRE bands for extended source assuming alpha=-1
# The result of hpXcalKcorr has to be multiplied by the effective beam area for alpha=-1
# NB: These are now read from SPIRE calibration tree
#k4E_PSW_new = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, BB=False, alpha=-1, \
#  ext=True, monoArea=beamMonoAreaPSW)[0] * beamAreaPipPSW
#k4E_PMW_new = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, BB=False, alpha=-1, \
#  ext=True, monoArea=beamMonoAreaPMW)[0] * beamAreaPipPMW
#k4E_PLW_new = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, BB=False, alpha=-1, \
#  ext=True, monoArea=beamMonoAreaPLW)[0] * beamAreaPipPLW

# Read K-correction factors for SPIRE bands for extended source assuming alpha=-1
k4E_PSW = fluxConv.meta['k4E_PSW'].double
k4E_PMW = fluxConv.meta['k4E_PMW'].double
k4E_PLW = fluxConv.meta['k4E_PLW'].double

kMonE_PSW = k4E_PSW / beamAreaPipPSW
kMonE_PMW = k4E_PMW / beamAreaPipPMW
kMonE_PLW = k4E_PLW / beamAreaPipPLW

# Calculate K-correction factors for Planck bands assuming alpha=-1
# although Planck uses extended sources, since its beam does not vary significantly
# with freqeuncy it is treated here as a point source
k4E_857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857,   BB=False, ext=False)[0]
k4E_545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545,   BB=False, ext=False)[0]

# Compute K-correction factors for SPIRE bands for extended source assuming alpha=-1
# NB: These are now read from SPIRE calibration tree
#k4P_PSW_new = hpXcalKcorr(spireRefFreq[0], freq, spireFiltPSW, BB=False, alpha=-1, ext=False)[0]
#k4P_PMW_new = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, BB=False, alpha=-1, ext=False)[0]
#k4P_PLW_new = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, BB=False, alpha=-1, ext=False)[0]

# Read K-correction factors for SPIRE bands for point source assuming alpha=-1
k4P_PSW = fluxConv.meta['k4P_PSW'].double
k4P_PMW = fluxConv.meta['k4P_PMW'].double
k4P_PLW = fluxConv.meta['k4P_PLW'].double

# Print Spire K4 factors for point source calibration for verification
print 'SPIRE k4 factors for point source:           %s, %s, %s' % (k4P_PSW,k4P_PMW,k4P_PLW)
#print 'CHECK: SPIRE k4 factors for point source:    %s, %s, %s' % (k4P_PSW_new,k4P_PMW_new,k4P_PLW_new)

# Print Spire K4 factors for extended source calibration for verification
print 'SPIRE k4 factors for extended source:        %s, %s, %s' % (k4E_PSW,k4E_PMW,k4E_PLW)
#print 'CHECK: SPIRE k4 factors for extended source: %s, %s, %s' % (k4E_PSW_new,k4E_PMW_new,k4E_PLW_new)

#-------------------------------------------------------------------------------
# Calculate and tabulate colour correction parameters for a range of temperatures
k545toPLW = Double1d()       # K-correction from HFI-545 to PLW
k857toPMW = Double1d()       # K-correction from HFI-857 to PMW
k857toPSW = Double1d()       # K-correction from HFI-857 to PSW
ratio545over857 = Double1d() # Ratio 545GHz to 857GHz filter

if fixBeta == True:
	#
	tvect = Double1d(range(ndiv)) * (Tmax - Tmin) / (ndiv-1) + Tmin
	for temp in tvect:
		#
		# K-corr for HFI-545 to monochromatic flux density
		k545 = hpXcalKcorr(hfiRefFreq[1], freq, hfiFilt545, \
		  BB=True, temp=temp, beta=beta0, ext=False)
		# K-corr for SPIRE PLW to monochromatic surface brightness
		kPLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, \
		  BB=True, temp=temp, beta=beta0, ext=True, monoArea=beamMonoAreaPLW)
		## convert to flux density using effective beam area
		##kPLW[0] = kPLW[0] * spireEffArea(freq, spireFiltOnlyPLW,
		#  beamMonoAreaPLW, BB=True, temp=temp, beta=beta0)
		## K-correction from HFI-545 to PLW		
		#k545toPLW.append(k545[0] / kPLW[0] * k4E_PLW / k4E_545 * kPLW[1] / k545[1])
		# K-correction from HFI-545 to PLW		
		k545toPLW.append(k545[0] / kPLW[0] * kMonE_PLW / k4E_545 * kPLW[1] / k545[1])
		#
		# K-corr for HFI-857 to monochromatic flux density
		k857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857, \
		  BB=True, temp=temp, beta=beta0, ext=False)
		# K-corr for SPIRE PMW to monochromatic surface brightness
		kPMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, \
		  BB=True, temp=temp, beta=beta0, ext=True, monoArea=beamMonoAreaPMW)
		# convert to flux density using effective beam area
		#kPMW[0] = kPMW[0] * spireEffArea(freq, spireFiltOnlyPMW,
		#  beamMonoAreaPMW, BB=True, temp=temp, beta=beta0)
		## K-correction from HFI-857 to PMW
		#k857toPMW.append(k857[0] / kPMW[0] * k4E_PMW / k4E_857 * kPMW[1] / k857[1])
		# K-correction from HFI-857 to PMW
		k857toPMW.append(k857[0] / kPMW[0] * kMonE_PMW / k4E_857 * kPMW[1] / k857[1])
		#
		# K-corr for SPIRE PMW to monochromatic surface brightness
		kPSW = hpXcalKcorr(spireRefFreq[0], freq,  spireFiltPSW, \
		  BB=True, temp=temp, beta=beta0, ext=True, monoArea=beamMonoAreaPSW)
		## convert to flux density using effective beam area
		#kPSW[0] = kPSW[0] * spireEffArea(freq, spireFiltOnlyPSW,
		#  beamMonoAreaPSW, BB=True, temp=temp, beta=beta0)
		## K-correction from HFI-857 to PSW
		#k857toPSW.append(k857[0] / kPSW[0] * k4E_PSW / k4E_857 * kPSW[1] / k857[1])
		# K-correction from HFI-857 to PSW
		k857toPSW.append(k857[0] / kPSW[0] * kMonE_PSW / k4E_857 * kPSW[1] / k857[1])
		#
		# Ratio of 545 and 845 GHz filters
		ratio545over857.append(k545[1] / k857[1] * k4E_545 / k4E_857 * k857[0] / k545[0])
	#
else:
	# Make beta vector
	bvect = Double1d(range(ndiv)) * (Bmax - Bmin) / (ndiv-1) + Bmin
	for beta in bvect:
		#
		# K-correction from HFI-545 to PLW
		k545 = hpXcalKcorr(hfiRefFreq[1],   freq, hfiFilt545, \
		  BB=True, temp=temp0, beta=beta, ext=False)
		kPLW = hpXcalKcorr(spireRefFreq[2], freq, spireFiltPLW, \
		  BB=True, temp=temp0, beta=beta, ext=True, monoArea=monoAreaPLW)
		# convert to flux density using effective beam area
		kPLW[0] = kPLW[0] * spireEffArea(freq, spireFiltOnlyPLW,
		  beamMonoAreaPLW, BB=True, temp=temp0, beta=beta)
		k545toPLW.append(k545[0] / kPLW[0] * k4E_PLW / k4E_545 * kPLW[1] / k545[1])
		#
		# K-correction from HFI-857 to PMW
		k857 = hpXcalKcorr(hfiRefFreq[0],   freq, hfiFilt857, \
		  BB=True, temp=temp0, beta=beta, ext=False)
		kPMW = hpXcalKcorr(spireRefFreq[1], freq, spireFiltPMW, \
		  BB=True, temp=temp0, beta=beta, ext=True, monoArea=monoAreaPMW)
		# convert to flux density using effective beam area
		kPMW[0] = kPMW[0] * spireEffArea(freq, spireFiltOnlyPMW,
		  beamMonoAreaPMW, BB=True, temp=temp0, beta=beta)
		k857toPMW.append(k857[0] / kPMW[0] * k4E_PMW / k4E_857 * kPMW[1] / k857[1])
		#
		# K-correction from HFI-857 to PSW
		kPSW = hpXcalKcorr(spireRefFreq[0], freq,  spireFiltPSW, \
		  BB=True, temp=temp0, beta=beta, ext=True, monoArea=monoAreaPSW)
		k857toPSW.append(k857[0] / kPSW[0] * k4E_PSW / k4E_857 * kPSW[1] / k857[1])
		# convert to flux density using effective beam area
		kPSW[0] = kPSW[0] * spireEffArea(freq, spireFiltOnlyPSW,
		  beamMonoAreaPSW, BB=True, temp=temp0, beta=beta)
		#
		# Ratio of 545 and 845 GHz filters
		ratio545over857.append(k545[1] / k857[1] * k4E_545 / k4E_857 * k857[0] / k545[0])
	#
#-------------------------------------------------------------------------------

if plot:
	# Plot colour correction factors
	if fixBeta == True:
		# First plot of parameters against dust temperature of assumed background
		p = PlotXY()
		p.addLayer(LayerXY(tvect,k545toPLW,name='k545toPLW'))
		p.addLayer(LayerXY(tvect,k857toPMW,name='k857toPMW'))
		p.addLayer(LayerXY(tvect,k857toPSW,name='k857toPSW'))
		p.addLayer(LayerXY(tvect,ratio545over857,name='ratio545over857'))
		p.xaxis.titleText = "T [K]"
		p.yaxis.titleText = "K-Factor"
		p.xaxis.type = Axis.LOG
		p.xaxis.range = [1.,100.]
		p.legend.visible = 1
	#
	else:
		# First plot of parameters against dust spectral index of assumed background
		p = PlotXY()
		p.addLayer(LayerXY(bvect,k545toPLW,name='k545toPLW'))
		p.addLayer(LayerXY(bvect,k857toPMW,name='k857toPMW'))
		p.addLayer(LayerXY(bvect,k857toPSW,name='k857toPSW'))
		p.addLayer(LayerXY(bvect,ratio545over857,name='ratio545over857'))
		p.xaxis.titleText = "Beta"
		p.yaxis.titleText = "K-Factor"
		p.legend.visible = 1
	#
	# Second plot of parameters against the ratio of the fluxes found in the 
	# all-sky maps of HFI
	p = PlotXY()
	p.addLayer(LayerXY(ratio545over857,k545toPLW,name='k545toPLW'))
	p.addLayer(LayerXY(ratio545over857,k857toPMW,name='k857toPMW'))
	p.addLayer(LayerXY(ratio545over857,k857toPSW,name='k857toPSW'))
	p.xaxis.titleText = "ratio 545/857"
	p.yaxis.titleText = "K-Factor"
	p.yaxis.range = [-0.1,4]
	p.legend.visible = 1

#-------------------------------------------------------------------------------
# Create calibration product with tabulated colour correction factors

# define the start and end dates for the product
# starting at beginning of PFM1
startDate = df.parse("2005.02.22/00:00:00/GMT")
endDate   = df.parse("2020.01.01/00:00:00/GMT")

photColorCorrHfi = herschel.spire.ia.dataset.PhotColorCorrHfi()

photColorCorrHfi.meta["creator"].value   = scriptVersionString
photColorCorrHfi.meta["modelName"].value = "FM"
photColorCorrHfi.meta["creationDate"].value = FineTime(java.util.Date())
photColorCorrHfi.meta["startDate"].value = FineTime(startDate)
photColorCorrHfi.meta["endDate"].value   = FineTime(endDate)
photColorCorrHfi.meta["fileOrigin"]  = herschel.ia.dataset.StringParameter(value="%s"%hfiRIMOFileName, description="Origin of the data")
photColorCorrHfi.setVersion(version)
photColorCorrHfi.setFormatVersion(formatVersion)

# Save standard point source correction factor into meta data
#photColorCorrHfi.meta["k4P_PSW"] = DoubleParameter(k4P_PSW, "PSW point source K-factor")
#photColorCorrHfi.meta["k4P_PMW"] = DoubleParameter(k4P_PMW, "PMW point source K-factor")
#photColorCorrHfi.meta["k4P_PLW"] = DoubleParameter(k4P_PLW, "PLW point source K-factor")
#photColorCorrHfi.meta["k4E_PSW"] = DoubleParameter(k4E_PSW, "PSW extended source K-factor")
#photColorCorrHfi.meta["k4E_PMW"] = DoubleParameter(k4E_PMW, "PMW extended source K-factor")
#photColorCorrHfi.meta["k4E_PLW"] = DoubleParameter(k4E_PLW, "PLW extended source K-factor")

# Save beta dependent colour correction tables into binary table
if fixBeta == True:
	photColorCorrHfi.meta["beta"] = DoubleParameter(beta0, "Modified black-body spectral index used")
        photColorCorrHfi.setTempVals(tvect)      # Temperature of modified BB
        #
# Save temperature dependent colour correction tables into binary table
else:
	photColorCorrHfi.meta['temperature'] = DoubleParameter(temp0, "Modified black-body temperature used")
	photColorCorrHfi.meta['temperature'].unit = Temperature.KELVIN
	#
	photColorCorrHfi.setTempVals(bvect)          # Spectral index of modified BB

photColorCorrHfi.meta["gamma"] = DoubleParameter(gamma, "Exponent describing FWHM dependence on frequency used")

# Save colour correction tables into binary table
photColorCorrHfi.setRatio545_857CorrVals(ratio545over857)  # Ratio 545/857 GHz filter
photColorCorrHfi.setK545toPLWCorrVals(k545toPLW)           # 545GHz to PLW K-factor
photColorCorrHfi.setK857toPMWCorrVals(k857toPMW)           # 857GHz to PMW K-factor
photColorCorrHfi.setK857toPSWCorrVals(k857toPSW)           # 857GHz to PSW K-factor


##########################################
# ******* write to FITS *******
#filename = java.io.File(r"%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_v%s.fits"%(directory, version))
if outputCalDirTree:
	filename = java.io.File(r"%s//Phot//SCalPhotColorCorrHfi//SCalPhotColorCorrHfi_v%s.fits"%(directory, version))
else:
	filename = java.io.File(r"%s//SCalPhotColorCorrHfi_v%s.fits"%(dataDir, version))
photColorCorrHfi.meta['fileName'] = herschel.ia.dataset.StringParameter(value=filename.name)
print
fitsWriter = FitsArchive()
fitsWriter.rules.append(herschel.spire.ia.util.MetaDataDictionary.getInstance().getFitsDictionary())
fitsWriter.save(filename.toString(), photColorCorrHfi)
print "written: %s"%filename.toString()
print

