#turns cal tables into tex tables for OM

cal=spireCal(pool='spire_cal_12_2')
spireBands=['PSW','PMW','PLW']
from herschel.ia.numeric.toolbox.util.MoreMath import modulo

#-------------------------------------------------------------------------------
##get Neptune beam areas
beamNepArc= {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamNeptunePswArc'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamNeptunePmwArc'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['beamNeptunePlwArc'].double}
##get Neptune alpha
alphaNep= {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['alphaNeptunePsw'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['alphaNeptunePmw'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['alphaNeptunePlw'].double}
##get effective frequencies
freqEff = {'PSW':cal.getPhot().getProduct("RadialCorrBeam").meta['freqEffPsw'].double, \
	'PMW':cal.getPhot().getProduct("RadialCorrBeam").meta['freqEffPmw'].double, \
	'PLW':cal.getPhot().getProduct("RadialCorrBeam").meta['freqEffPlw'].double}

#-------------------------------------------------------------------------------
##get beam correction factors
kBeam = cal.getPhot().getProduct("ColorCorrBeam")

##get pipeline beam areas
pipBeamArc = {'PSW':cal.getPhot().getProduct("ColorCorrBeam").meta['beamPswArc'].double, \
	'PMW':cal.getPhot().getProduct("ColorCorrBeam").meta['beamPmwArc'].double, \
	'PLW':cal.getPhot().getProduct("ColorCorrBeam").meta['beamPlwArc'].double}

pipBeamSr = {'PSW':cal.getPhot().getProduct("ColorCorrBeam").meta['beamPswSr'].double, \
	'PMW':cal.getPhot().getProduct("ColorCorrBeam").meta['beamPmwSr'].double, \
	'PLW':cal.getPhot().getProduct("ColorCorrBeam").meta['beamPlwSr'].double}

#calculate effective beam areas
effBeamArc = kBeam.copy()
for n in kBeam.getSets():
	print n
	for band in spireBands:
		print band
		effBeamArc[n][band].data = pipBeamArc[band]/kBeam[n][band].data

#-------------------------------------------------------------------------------
##get kColP, kColE tables
kPsrc=cal.getPhot().getProduct("ColorCorrKList")[1]
kExtd=cal.getPhot().getProduct("ColorCorrKList")[0]

#-------------------------------------------------------------------------------
##get k4 parameters
k4P={'PSW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4P_PSW'].double,\
	'PMW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4P_PMW'].double, \
	'PLW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4P_PLW'].double}
k4E={'PSW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4E_PSW'].double,\
	'PMW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4E_PMW'].double, \
	'PLW':cal.getPhot().getProduct("FluxConvList")[0].meta['k4E_PLW'].double}

kMonE={}
kPtoE={}
k4E4P={}
for band in spireBands:
	kMonE[band] = (k4E[band] / pipBeamSr[band]) / 1.e6
	kPtoE[band] = (kMonE[band] / k4P[band])
	k4E4P[band] = k4E[band] / k4P[band]

#-------------------------------------------------------------------------------
##get aperture corrections
apCorr_incBG = cal.getPhot().getProduct('ColorCorrApertureList')[0]
apCorr_noBG = cal.getPhot().getProduct('ColorCorrApertureList')[1]

################################################################################
## Print relevant OM tables
################################################################################

print '\n----- Table 5.2 (Basic 2-D Gaussian parameters), last 2 rows:'
print '\\begin{tabular}{l|ccc}'
print '\\hline\\hline'
print 'Band & PSW & PMW & PLW\\\\'
print '\hline'
print '$\\alpha_\\mathrm{Nep}$ & %.2f & %.2f & %.2f \\\\'%(alphaNep['PSW'],alphaNep['PMW'],alphaNep['PLW'])
print 'Major$\\times$Minor\\- FWHM (arcsec) & 18.3$\\times$17.0 & 24.7$\\times$23.2 & 37.0$\\times$33.4 \\\\'
print 'Geometric mean FWHM ($\\theta_\\mathrm{Nep}$, arcsec) & 17.6 & 23.9 & 35.2 \\\\'
print 'Ellip\\-ticity (\\%%)& 8.1 & 6.6 & 10.9 \\\\'
print 'Measured beam solid angle ($\\Omega_\\mathrm{Nep}$, arcsec$^2$) & 450 & 795 & 1665 \\\\'
print 'Measured beam solid angle ($\\Omega_\\mathrm{Nep}$, arcsec$^2$) & %.0f & %.0f & %.0f \\\\'%\
	(beamNepArc['PSW'],beamNepArc['PMW'],beamNepArc['PLW'])
print 'Isophotal frequency$^a$ ($\\nu_\\mathrm{eff}$, GHz) & %.2f & %.2f & %.2f'%\
	(freqEff['PSW'],freqEff['PMW'],freqEff['PLW'])
print '\\end{tabular}'

#-------------------------------------------------------------------------------
print '\n----- Table 5.3 (Effective beam solid angle, alpha):'
print '\\begin{tabular}{r|ccc|ccc}'
print '\\hline\\hline'
print '& \\multicolumn{3}{c|}{Effective beam solid angle} & \\multicolumn{3}{c}{Beam correction factor} \\\\'
print '& \\multicolumn{3}{c|}{$\\Omega_\\mathrm{eff}$ (arcsec$^2$)} & \\multicolumn{3}{c}{$K_\\mathrm{Beam}$} \\\\'
print '$\\alpha_S$  & PSW & PMW & PLW & PSW & PMW & PLW \\\\'
print '\\hline'
kB=kBeam['alpha']
eB=effBeamArc['alpha']
for a in range(len(eB['alpha'].data)):
	if eB['alpha'].data[a]==-1.:
		print '{\\bf %.1f }& {\\bf %.2f} & {\\bf %.2f} & {\\bf %.2f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} \\\\'%\
			(eB['alpha'].data[a],\
			eB['PSW'].data[a],eB['PMW'].data[a],eB['PLW'].data[a],\
			kB['PSW'].data[a],kB['PMW'].data[a],kB['PLW'].data[a])
	else:
		print '%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\'%\
			(eB['alpha'].data[a],\
			eB['PSW'].data[a],eB['PMW'].data[a],eB['PLW'].data[a],\
			kB['PSW'].data[a],kB['PMW'].data[a],kB['PLW'].data[a])
print '\\hline'
print '\\end{tabular}'

#-------------------------------------------------------------------------------
print '\n----- Table 5.4 (Effective beam solid angle, beta=1.5,2.0):'
print '\\begin{tabular}{r|ccc|ccc}'
print '\\hline\\hline'
print '& \\multicolumn{3}{c|}{Effective beam solid angle} & \\multicolumn{3}{c}{Beam correction factor} \\\\'
print '& \\multicolumn{3}{c|}{$\\Omega_\\mathrm{eff}(T,\\beta)$ (arcsec$^2$)} & \\multicolumn{3}{c}{$K_\\mathrm{Beam}(T,\\beta)$} \\\\'
print 'Temp (K) & PSW & PMW & PLW & PSW & PMW & PLW \\\\'
kB=kBeam['beta_1_50']
eB=effBeamArc['beta_1_50']
print '\\hline'
print ' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\'
print '\\hline'
for t in range(len(eB['Temperature'].data)):
	if (eB['Temperature'].data[t] < 10) or (eB['Temperature'].data[t] <= 40 and modulo(eB['Temperature'].data[t],5)==0.):
		print '%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\'%\
			(eB['Temperature'].data[t],\
			eB['PSW'].data[t],eB['PMW'].data[t],eB['PLW'].data[t],\
			kB['PSW'].data[t],kB['PMW'].data[t],kB['PLW'].data[t])
kB=kBeam['beta_2_00']
eB=effBeamArc['beta_2_00']
print '\\hline'
print ' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\'
print '\\hline'
for t in range(len(eB['Temperature'].data)):
	if (eB['Temperature'].data[t] < 10) or (eB['Temperature'].data[t] <= 40 and modulo(eB['Temperature'].data[t],5)==0.):
		print '%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\'%\
			(eB['Temperature'].data[t],\
			eB['PSW'].data[t],eB['PMW'].data[t],eB['PLW'].data[t],\
			kB['PSW'].data[t],kB['PMW'].data[t],kB['PLW'].data[t])
print '\\hline'
print '\\end{tabular}'

#-------------------------------------------------------------------------------
print '\n----- Table 5.5 (KColP, KColE with alpha):'
print '\\begin{tabular}{c|lll|lll|}'
print '\\hline\\hline'
print '& \\multicolumn{3}{c|}{Point Source ($K_\\mathrm{ColP}$)} & \\multicolumn{3}{c|}{Extended source ($K_\\mathrm{ColE}$)} \\\\'
print '$\\alpha_S$  & PSW & PMW & PLW & PSW & PMW & PLW \\\\'
print '\\hline'
kP=kPsrc['alpha']
kE=kExtd['alpha']
for a in range(len(kP['alpha'].data)):
	if kP['alpha'].data[a]==2.:
		#bold print
		print '{\\bf %.1f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} & {\\bf %.4f} \\\\'%\
			(kP['alpha'].data[a],\
			kP['PSW'].data[a],kP['PMW'].data[a],kP['PLW'].data[a],\
			kE['PSW'].data[a],kE['PMW'].data[a],kE['PLW'].data[a])
	else:
		print '%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\'%\
			(kP['alpha'].data[a],\
			kP['PSW'].data[a],kP['PMW'].data[a],kP['PLW'].data[a],\
			kE['PSW'].data[a],kE['PMW'].data[a],kE['PLW'].data[a])
print '\\end{tabular}'

#-------------------------------------------------------------------------------
print '\n----- Table 5.6 (KColP, kColE with beta=1.5,2.0):'
print '\\begin{tabular}{r|lll|lll}'
print '\\hline\\hline'
print '& \\multicolumn{3}{c|}{Point source ($K_\\mathrm{ColP}$)} &'
print '\\multicolumn{3}{c}{Extended source ($K_\\mathrm{ColE}$)} \\\\'
print '$T$ & PSW & PMW & PLW & PSW & PMW & PLW \\\\'
kP=kPsrc['beta_1_50']
kE=kExtd['beta_1_50']
print '\\hline'
print ' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\'
print '\\hline'
for t in range(len(eB['Temperature'].data)):
	if (kP['Temperature'].data[t] < 10) or (kP['Temperature'].data[t] <= 40 and modulo(kP['Temperature'].data[t],5)==0.):
		print '%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\'%\
			(kP['Temperature'].data[t],\
			kP['PSW'].data[t],kP['PMW'].data[t],kP['PLW'].data[t],\
			kE['PSW'].data[t],kE['PMW'].data[t],kE['PLW'].data[t])
kP=kPsrc['beta_2_00']
kE=kExtd['beta_2_00']
print '\\hline'
print ' & \multicolumn{6}{c}{$\\beta=1.5$}\\\\'
print '\\hline'
for t in range(len(eB['Temperature'].data)):
	if (kP['Temperature'].data[t] < 10) or (kP['Temperature'].data[t] <= 40 and modulo(kP['Temperature'].data[t],5)==0.):
		print '%.1f & %.2f & %.2f & %.2f & %.4f & %.4f & %.4f \\\\'%\
			(kP['Temperature'].data[t],\
			kP['PSW'].data[t],kP['PMW'].data[t],kP['PLW'].data[t],\
			kE['PSW'].data[t],kE['PMW'].data[t],kE['PLW'].data[t])
print '\\hline'
print '\\end{tabular}'

#-------------------------------------------------------------------------------
print '\n----- Table 5.7 (k4 parameters etc.):'
print '\\begin{tabular}{l|ccc}'
print '\\hline\\hline'
print '& PSW & PMW & PLW \\\\'
print '\\hline'
print '$K_\\mathrm{4P}$ & %.4f & %.4f & %.4f \\\\'%(k4P['PSW'],k4P['PMW'],k4P['PLW'])
print '$K_\\mathrm{MonE}$ (MJy/sr per Jy/beam) & %.3f & %.3f & %.3f \\\\'%(kMonE['PSW'],kMonE['PMW'],kMonE['PLW'])
print '$K_\\mathrm{PtoE}$ (MJy/sr per Jy/beam) & %.3f & %.3f & %.3f \\\\'%(kPtoE['PSW'],kPtoE['PMW'],kPtoE['PLW'])
print '$\\Omega_\\mathrm{pip}$ (arcsec$^2$) & %.2f & %.2f & %.2f \\\\'%(pipBeamArc['PSW'],pipBeamArc['PMW'],pipBeamArc['PLW'])
print '\\hline'
print '$K_\\mathrm{4E}$ & %.4f & %.4f & %.4f \\\\'%(k4E['PSW'],k4E['PMW'],k4E['PLW'])
print '$K_\\mathrm{4E}/K_\\mathrm{4P}$ & %.4f & %.4f & %.4f '%(k4E4P['PSW'],k4E4P['PMW'],k4E4P['PLW'])
print '\\end{tabular}'

#-------------------------------------------------------------------------------
print '\n----- Table 5.8 (Aperture correction):'
print'\\begin{tabular}{r|ccc|ccc}'
print '\\hline\\hline'
print '& \\multicolumn{3}{c|}{Background included} & \\multicolumn{3}{c}{No Background} \\\\'
print '$\\alpha$ & PSW & PMW & PLW & PSW & PMW & PLW \\\\'
print '\\hline'
aI=apCorr_incBG['alpha']
aN=apCorr_noBG['alpha']
for a in range(len(aI['alpha'].data)):
	print '%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\'%\
	(aI['alpha'].data[a],\
	aI['PSW'].data[a],aI['PMW'].data[a],aI['PLW'].data[a], \
	aN['PSW'].data[a],aN['PMW'].data[a],aN['PLW'].data[a])
print '\\end{tabular}'
