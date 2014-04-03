#-------------------------------------------------------------------------------
#===============================================================================
#=====                           DEFINE FUNCTIONS                          =====
#===============================================================================
#-------------------------------------------------------------------------------

# Set up functions to calculate beam profile, effective frequency and effective area
# Function list:
#   * spireMonoBeam: Calculate monochromatic beam profile and area at a given frequency
#   * spireMonoAreas: Calculate monochromatic beam areas over a range of frequencies
#   * spireMonoAreasSimple: Calculate monochromatic beam areas over a range of frequencies
#        using simple beam model [beam area ~ (freq/effFreq)**2gamma]
#   * spireEffArea: Calculate the effective beam area for a given spectrum
#   * spireFindEffFreq: Calculate the effective frequency for SPIRE
#   * spireEffBeam: Calculate the effective beam profile, area and beam map for SPIRE
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
	               beams areas should be Beam Radial Profiles calculated
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
def spireEffArea(freq, transm, monoArea, verbose=False, BB=False, UD=False, \
                 alpha=-1.0, temp=20.0, beta=1.8, userSpectrum=None):
	"""
	========================================================================
	Calculate the effective beam area for a source of a given spectrum

	Inputs:
	  freq:       (array float) frequency vector corresponding to RSRF values [Hz]
	  transm:     (array float) relative spectral response (RSRF) corresponding to freq
	                Note that this should *not* include the aperture efficiency
	  monoArea:   (array float) monochromatic beam solid angle corresponding
	                to frequencies in freq
          verbose:    (boolean) Set to print additional information.
                        OPTIONAL. Default=False
	  BB:         (boolean) set to use modified black body spectrum with
                        temperature temp and emissivity index beta
	                OPTIONAL. Default=False
	  UD:         (boolean) set to use user-defined spectrum given by
                        userSpec
	                OPTIONAL. Default=False
	  alpha:      (float) Exponent of power-law sky background model
                        (if BB==False and UD==False)
	                OPTIONAL. Default=-1
	  temp:       (float) Dust/sky temperature (if BB=True)
	                OPTIONAL. Default=20.0
	  beta:       (float) Dust/sky spectral index (if BB=True)
	                OPTIONAL. Default=1.8
	  userSpec:   (array float) User-defined spectrum corresponding to freq.
                        (if UD=True)
	                OPTIONAL. Default=None
	Outputs:     
	            (float) Beam area for given spectrum, in same units as monoArea

	Calculation:
	  Calculates the source spectrum (either modifies black body, power law
            or user-defined)
          Multiplies the monochromatic beam area by RSRF and source spectrum
	  Integrates over frequency
	  Normalises by integral over frequency of RSRF and source spectrum

	Dependencies:
	  herschel.ia.numeric.toolbox.interp.LinearInterpolator
	  herschel.ia.numeric.toolbox.integr.TrapezoidalIntegrator

	2013/12/19  C. North  initial version

	"""	

        # check not both BB and UD are true
        assert (BB and UD)==False,'ERROR: both BB and UD are set'
        #check monoArea and transm are of right length
        assert len(monoArea)==len(freq),\
                'ERROR: monoArea myst be of same length as freq (%d,%d)'%\
                (len(monoArea),len(freq))
        assert len(transm)==len(freq),\
                'ERROR: transm myst be of same length as freq (%d,%d)'%\
                (len(transm),len(freq))
	#
	# Calculate sky background model
	#
	if BB == 1:
                # Modified black body
		if verbose:
                        print 'Using greybody with T=%.2f, beta=%.3f'%(temp,beta)
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
        elif UD == 1:
                # User-defined spectrum
                # Check userSpec is of same length
                assert len(userSpec)==len(freq),\
                        'ERROR: userSpec must be of same length as freq (%d,%d)'%(len(userSpec),len(freq))
		if verbose:
                        print 'Using user-defined spectrum'
                fSky=userSpec
                
	else:
                # Power-Law spectrum
		if verbose:
                        print 'Using power law spectrum with spectral index %.2f'%(alpha)
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
# Calculate the effective frequency for SPIRE
def spireFindEffFreq(freq, rsrf, beamProfs, effFreqInit, gamma,
  areaNep, alphaNep, array, freqFact=500,
  initRange=0.01, reqPrec=1.e-6, maxIter=5, verbose=False):
	"""
	========================================================================
	Derive effective frequency for spire bands. This is the frequency at
	which the monochromatic beam area is equal to the area as measured on
	Neptune.

	Inputs:
	  freq:        (array float) frequency vector corresponding to RSRF values [Hz]
	  rsrf:        (array float) relative spectral response (RSRF) corresponding to freq
	                 Note that this should *not* include the aperture efficiency
	  beamProfs:   (dataset) PhotRadialCorrBeam object from calibration tree
	  effFreqInit: (float) Initial estimate of effective frequency [Hz]
	  gamma:       (float) Exponent of powerlaw describing FWHM dependence
	                 on frequency
	  areaNep:     (float) Solid angle of beam measured on Neptune [sr]
	  alphaNep:    (float) Spectral index of Neptune frequency spectrum
	  array:       (string) spire array ('Psw'|'Pmw'|'Plw')
	  freqFact:    (int) Factor by which to reduce size of freq for calculations.
	                 Only applicable if SimpleBeam=False
	                 OPTIONAL. Default=500.
	  initRange:   (float) Fractional intial range to use in calculations
	                 OPTIONAL. Default=0.01
	  reqPrec:     (float) Required relative precision for convergence.
	                 OPTIONAL. Default=1.e-6
          maxIter:     (int) Maximum interations to try
	                 OPTIONAL. Default=5
	  verbose:     (boolean) set to print more detailed info
                         OPTIONAL. Default=False

	Outputs:
	               (float) Effective frequency [Hz]

	Calculation:
	  Uses effFreqInit +/- initRange to calculate monochromatic Areas
	  Calculates effective beam area for Neptune spectrum (alphaNep)
 	  Compares with measured Neptune beam area (areaNep)
	  Adjusts effFreqInit and iterates until maximum iterations (maxIter) or
	  required precicion (reqPrec) is reached

	Dependencies:
	  spireMonoAreasSimple
	  spireMonoAreas
	  spireEffArea

	2014/01/07  C. North  initial version
	"""
		
	#calculate initial estimates of effective Frequency
	#parameterised be offset from original estimate
	relEff=Double1d([1.-initRange,1.+initRange])
	effFreqs=effFreqInit*relEff
	if verbose:
		print 'Calculating effective frequency for %s'%array
		print '  Initial %s Effective Frequencies: [%.2f : %.2f] GHz'%(array,effFreqs[0]/1.e9,effFreqs[1]/1.e9)
	#calculate effective beam area for initial estimates of effFreq
        beamMonoArea0=spireMonoAreas(freq,beamProfs,effFreqs[0],gamma,array,freqFact=freqFact)
        beamMonoArea1=spireMonoAreas(freq,beamProfs,effFreqs[1],gamma,array,freqFact=freqFact)
	beamMonoDiff=beamMonoArea1-beamMonoArea0
	effAreas=Double1d(2)
	effAreas[0]=spireEffArea(freq,rsrf, beamMonoArea0, BB=False, alpha=alphaNep)
	effAreas[1]=spireEffArea(freq,rsrf, beamMonoArea1, BB=False, alpha=alphaNep)

	iter=0
	done=False
	while ((done==False) and (iter <= maxIter)):
		iter=iter+1
		#difference from measured beam area
		diffAreas=effAreas-areaNep
		relAreas=diffAreas/areaNep
		#calculate new esitmate of rel
		grad=(diffAreas[1]-diffAreas[0])/(relEff[1]-relEff[0])
		relEffNew=relEff[1] - diffAreas[1]/grad

		#move values in arrays
		relEff[0]=relEff[1]
		effAreas[0]=effAreas[1]

		#calculate new effective beam area
		relEff[1]=relEffNew
		effFreqs=effFreqInit*relEff

                beamMonoNew=spireMonoAreas(freq,beamProfs,effFreqs[1],gamma,array,freqFact=freqFact)
		effAreas[1]=spireEffArea(freq,rsrf, beamMonoNew, BB=False, alpha=alphaNep)

		diffAreas=effAreas-areaNep
		relAreas=diffAreas/areaNep
		if verbose:
			print '    iter %d: %.4f [effFreq %.2f GHz], Area=%.2f, RelDiff=%.4g'%(iter,relEff[1],effFreqs[1]/1.e9,effAreas[1]/arcsec2Sr,relAreas[1])
		if (Math.abs(relAreas[1]) < reqPrec):
			done=True
	if ((iter > maxIter) and (done==False)):
		print "  Warning: maximum iterations [%d] exceeded without conversion [%g]"%(maxIter,reqPrec)

	if verbose:
		print '  Final %s effFreq: %.4f'%(array,effFreqs[1]/1.e9)
		print '  Resulting %s Neptune area: %.2f [rel Diff: %.3g]'%(array,effAreas[1]/arcsec2Sr,relAreas[1])
	return(effFreqs[1])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Calculate the effective beam profile and area for SPIRE
def spireEffBeam(freq, transm, beamProfs, effFreq, gamma, array,
  BB=False,UD=False,temp=20.0,beta=1.8,alpha=-1.0,userSpec=None,verbose=False):

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
	  BB:         (boolean) set to use modified black body spectrum with
                        temperature temp and emissivity index beta
	                OPTIONAL. Default=False
	  UD:         (boolean) set to use user-defined spectrum given by
                        userSpec
	                OPTIONAL. Default=False
	  alpha:      (float) Exponent of power-law sky background model
                        (if BB==False and UD==False)
	                OPTIONAL. Default=-1
	  temp:       (float) Dust/sky temperature (if BB=True)
	                OPTIONAL. Default=20.0
	  beta:       (float) Dust/sky spectral index (if BB=True)
	                OPTIONAL. Default=1.8
	  userSpec:   (array float) User-defined spectrum corresponding to freq.
                        (if UD=True)
	                OPTIONAL. Default=None
          verbose:    (boolean) Set to print additional information.
                        OPTIONAL. Default=False.

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

        # check not both BB and UD are true
        assert (BB and UD)==False,'ERROR: both BB and UD are set'
        #check monoArea and transm are of right length
        assert len(monoArea)==len(freq),\
                'ERROR: monoArea myst be of same length as freq (%d,%d)'%\
                (len(monoArea),len(freq))
        assert len(transm)==len(freq),\
                'ERROR: transm myst be of same length as freq (%d,%d)'%\
                (len(transm),len(freq))
	#
	# Calculate sky background model
	#
	if BB == 1:
                # Modified black body
		if verbose:
                        print 'Using greybody with T=%.2f, beta=%.3f'%(temp,beta)
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
	#
        elif UD == 1:
                # User-defined spectrum
                # Check userSpec is of same length
                assert len(userSpec)==len(freq),\
                        'ERROR: userSpec must be of same length as freq (%d,%d)'%(len(userSpec),len(freq))
		if verbose:
                        print 'Using user-defined spectrum'
                fSky=userSpec
                
	else:
                # Power-Law spectrum
		if verbose:
                        print 'Using power law spectrum with spectral index %.2f'%(alpha)
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
# Calculate K-correction parameters for given spectrum & source type
def hpXcalKcorr(freq0, freq, transm, BB=True, UD=True, temp=20.0, beta=1.8, alpha=-1.0,
                userSpec=None, ext=False, monoArea=None, verbose=False):
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
	  BB:         (boolean) set to use modified black body spectrum with
                        temperature temp and emissivity index beta
	                OPTIONAL. Default=False
	  UD:        (boolean) set to use user-defined spectrum given by
                       userSpec
	               OPTIONAL. Default=False
	  alpha:     (float) Exponent of power-law sky background model
                       (if BB==False and UD==False)
	               OPTIONAL. Default=-1
	  temp:      (float) Dust/sky temperature (if BB=True)
	               OPTIONAL. Default=20.0
	  beta:      (float) Dust/sky spectral index (if BB=True)
	               OPTIONAL. Default=1.8
	  userSpec:  (array float) User-defined spectrum corresponding to freq.
                       (if UD=True)
	               OPTIONAL. Default=None
	  ext:       (boolean) calculating for extended source
	                OPTIONAL. Default=False
	  monoArea:  (array float) Monochromatic Beam solid angle [Sr] corresponding
	                to freq.
	                OPTIONAL. Only required if ext=True
          verbose:   (boolean) Set to print more information.
                        OPTIONAL. Default=False.

	Outputs:
	 (double array) [0]: K-correction factor
	                [1]: Sky emission at reference fequency (fSky0)
	
	Calculation:
	  Depending on the state of the input parameter BB and UD, uses spectrum of a
	  Planck function multiplied by frequency to the power of beta, a user-defined
          spectrum, or a power-law spectrum with spectral index alpha. The spectrum
          is calculated for all values in the vector freq. In addition the same value
          is calculated at the discrete frequency freq0. Then the product of this 
          value and the integral of the RSRF over all frequencies, divided by the 
          integral over all products of frequency and RSRF is calculated.

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
        # check not both BB and UD are true
        assert (BB and UD)==False,'ERROR: both BB and UD are set'
        #check monoArea and transm are of right length
        assert len(monoArea)==len(freq),\
                'ERROR: monoArea myst be of same length as freq (%d,%d)'%\
                (len(monoArea),len(freq))
        assert len(transm)==len(freq),\
                'ERROR: transm myst be of same length as freq (%d,%d)'%\
                (len(transm),len(freq))

	#
	# Calculate sky background model
	#
	# 1) As a modified Blackbody
	if BB == 1:
		if verbose:
                        print 'Using greybody with T=%.2f, beta=%.3f'%(temp,beta)
		fSky  = 2*h * freq**3 / c**2 / (EXP(h*freq/k/temp) - 1.) * freq**beta
		fSky0 = 2*h * freq0**3 / c**2 / (EXP(h*freq0/k/temp) - 1.) * freq0**beta
	#
        elif UD == 1:
                # User-defined spectrum
                # Check userSpec is of same length
                assert len(userSpec)==len(freq),\
                        'ERROR: userSpec must be of same length as freq (%d,%d)'%(len(userSpec),len(freq))
		if verbose:
                        print 'Using user-defined spectrum'
                fSky=userSpec
                fsKy0=CubicSplineInterplator(freq,userSpec)(freq0)
	else:
                # Power-Law
		if verbose:
                        print 'Using power law spectrum with spectral index %.2f'%(alpha)
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

