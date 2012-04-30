#
# Mercury mapping project binning the wetlands data
# a m!ke@usgs joint
#
# mnfienen@usgs.gov
#
#

import rpy2                     # bring in the general rpy2 module
import rpy2.rinterface as rinterface # to convert R objects into python
import rpy2.robjects as rob     # bring in the robjects module as r since it is used the most
import rpy2.rlike.container as rlc # used to make data.frame objects
import numpy as np              # no work is complete without NUMPY!
import rpy2.robjects.numpy2ri   # this makes numpy -> rpy conversions possible
r = rob.r
import scat_plotter as sp       # scatter plotting routine
import os as os
import sys
import pickle


from predplot_especiale import mypredplot


rpy2.robjects.numpy2ri.activate()

######
# function to convert a list to a 2d float array
# 1) converts list to a 2-d array
# 2) if requested by logtrans flag, performs log transformation
######
def lst2array_float(inlst, logtrans):
    inlst = np.array(inlst).astype(float)
    inlst = np.atleast_2d(inlst)
    if logtrans:
        inlst = np.log10(inlst)
    inlst[inlst==-999]=np.NaN
    return inlst


    
    

########################################################################
                              # MAIN # 
########################################################################
# number of optimized components using one_Standard_error criterion

opt_comp = 3

###########
# PREDICTIONS FLAG - IF TRUE, INVOKE PREDICTION CAPABILITIES AT THE END!
PredFlag = False
pred_infile = 'mehg_prediction_input.dat'

###########
# flag to determine how replicates used.  True means replications are averaged
reps_avg_flag = True

# Calibration Data File names
if reps_avg_flag:
    infile = 'hg_regressionAvgReps.dat'
else:
    infile = 'hg_regression_studies.dat'

# set up control flags for transformations, etc.
HG_RAT_FLAG = False    # ratio flag for methyl/total mercury True means yes ratio
                      # false means Methyl only
                      
# default path rootname
newpath = ''
# specify the number of variables to use
nvars           = 3 # CURRENTLY ONLY FOR OUTPUT CONTROL --> changed to 4 if use wetlands

# wetlflag --> flag for whether to use flag for wetlands
wetlflag        = True
if wetlflag:
    nvars = 4

# reducing dataset due to drainage size bounds.  
# N.B. this used to be optional (through a flag) but now it's required.  To disable, just make the
# upper bound greater than the largest in the dataset!
upper = np.atleast_2d(np.array([250000])).T #these are the upper bounds to consider
lower = np.zeros_like(upper)
db = np.hstack((lower,upper))

# the following flags are for log10 transformations.  True means TRANSFORM
HGlogFlag       = True
so4logFlag      = True
toclogFlag      = True
wetlandlogFlag  = True
#
############

############
# Individual study names
# 0-LNWR_FINAL_DATA
# 1-NAWQA_FL_WI_STREAMWATER_OR
# 2-WLS_HG_MEHG_ANCILLARY
# 3-NAWQA_NY_SC_H20SELECT_20090206
# 4-HG_MINNESOTANWIS_20090121
# 5-NPS_LAKES_WOUT_MILLS
# 6-NAWQA_HGNATLSYN_BCS_ENV_DATA
# 7-INDIANAHGSTREAMS2007_2008
allstudies = ['LNWR_FINAL_DATA','NAWQA_FL_WI_STREAMWATER_OR','WLS_HG_MEHG_ANCILLARY', \
              'NAWQA_NY_SC_H20SELECT_20090206', 'HG_MINNESOTANWIS_20090121', \
              'NPS_LAKES_WOUT_MILLS', 'NAWQA_HGNATLSYN_BCS_ENV_DATA', \
              'INDIANAHGSTREAMS2007_2008']
# Select current study to consider as an index from above.
# For all, cstudy = 99999
cstudy=99999
if (cstudy < 9999):
    currstud = allstudies[cstudy]
#
############


for drainbounds in db:

	# set up empty bins for the data
	siteid      = []
	toc         = []
	so4         = []
	pH          = []
	methHg      = []
	totHg       = []
	wetland     = []
	cumdrainage = []
	if (cstudy < 9999):
		study   = []
	
	
	# read in the data table  and instantly convert to an array
	for line in open(infile,'r'):
		tmpline = line.strip().split('~')
		siteid.append(tmpline[0])
		toc.append(tmpline[3])
		so4.append(tmpline[4])
		pH.append(tmpline[5])
		methHg.append(tmpline[6])
		totHg.append(tmpline[7])
		wetland.append(tmpline[17])
		cumdrainage.append(tmpline[8])
		if (cstudy < 9999):
			study.append(tmpline[18])
		
	#make arrays of the float variables and logtrans as appropriate
	toc     = lst2array_float(toc,toclogFlag)
	so4     = lst2array_float(so4,so4logFlag)
	pH      = lst2array_float(pH,False)
	methHg  = lst2array_float(methHg,HGlogFlag)
	totHg   = lst2array_float(totHg,HGlogFlag)
	wetland = lst2array_float(wetland,wetlandlogFlag)
	cumdrainage = lst2array_float(cumdrainage,False)
	if (cstudy < 9999):
		study = np.asarray(study).astype(str)
	
	# define the dependent variable as either methyl mercury or the fraction
	#    of methyl mercury to total
	if HG_RAT_FLAG:
			# account for log-space division
		if HGlogFlag:
			HG_reg = methHg-totHg
		else:
			HG_reg = methHg/totHg
		hglabel = 'Methyl Hg/Total Hg'
	else:
		HG_reg = methHg        
		hglabel = 'Methyl Hg'
		
	# kludge fix for log of zero
	if wetlandlogFlag:
		wetland[np.isinf(wetland)] = -3.5
	if toclogFlag:
		toc[np.isinf(toc)] = -2
	
	# form the RHS matrix and dependent vector
	if wetlflag:
		RHS = np.hstack((so4.T,toc.T,pH.T,wetland.T))
	else:
		RHS = np.hstack((so4.T,toc.T,pH.T))
	HG_reg = HG_reg.T
	
	#
	# if a specific study is indicated, purge all else from the data
	#
	if (cstudy < 9999):
		inds_to_keep = np.nonzero(study==currstud)
		HG_reg = HG_reg[inds_to_keep]
		RHS = RHS[inds_to_keep,:].squeeze()

	inds_to_keep = np.nonzero(np.logical_and(cumdrainage>=drainbounds[0],cumdrainage<=drainbounds[1])==True)[1]
	HG_reg = HG_reg[inds_to_keep]
	RHS = RHS[inds_to_keep,:].squeeze()
	newpath = 'cd_lb' + str(int(drainbounds[0])) + '_ub' + str(int(drainbounds[1])) 
	boundname = newpath
	if (not wetlflag):
		newpath += '_noWL_'
	if (os.path.exists(newpath)==False):
		os.mkdir(newpath)
	newpath += '/'
		
	xname = ['SULFATE','TOC','pH','Wetland %']
	
	for i in xrange(RHS.shape[1]):
		xvals = RHS[:,i]
		yvals = HG_reg
		title = 'Scatter ' + xname[i]
		if (cstudy < 9999):
			outfile = xname[i] + '_' + currstud + '_scatter.png'
		else:
			outfile = xname[i] + '_scatter.png'
		outfile = newpath  + outfile
		showme = False
		sp.scat_plot(i,xvals,yvals,title,xname[i],outfile,HG_RAT_FLAG,showme)
	
	if (cstudy < 9999):
		outHGfile = 'methyl_to_total_HG'+ '_' + currstud + '.png'   
	else:
		outHGfile = 'methyl_to_total_HG.png'
	outHGfile = newpath  + outHGfile
	sp.scat_plot(10,np.arange(1,len(HG_reg)+1),HG_reg,hglabel,'meas. number',outHGfile,HG_RAT_FLAG,showme)
	
	
	HGr     = rob.conversion.py2ri(HG_reg)
	RHSr    = rob.conversion.py2ri(RHS)
	#
	CalData = rlc.TaggedList([HGr,RHSr],tags=('hg','rhs'))
	CalData = rob.DataFrame(CalData)
	
	r('''library(pls)''')
	#rob.globalEnv["HGr"] = HGr
	#rob.globalEnv["RHSr"] = RHSr
	rob.globalenv["CalData"] = CalData
	
	
	# perform the PLS regression
	if wetlflag:
		HGresults = r.plsr(r("hg ~ rhs.1 + rhs.2 + rhs.3 + rhs.4"),data=CalData,validation="LOO")
	else:
		HGresults = r.plsr(r("hg ~ rhs.1 + rhs.2 + rhs.3"),data=CalData,validation="LOO")

	# Save down (pickle) the results for optimization of components later
	# first extract the predictions from the cross-validation
	caseID = os.getcwd().split('/')[-1]

	nam1 = list(HGresults.names)
	cv_d1 = HGresults[nam1.index('validation')]
	nam2 = list(cv_d1.names)
	cv_data = np.array(cv_d1[nam2.index('pred')]).squeeze()
	ofp = open(caseID + '_' + 'cv_data.pkl','wb')
	pickle.dump(cv_data,ofp)
	ofp.close()
	# then save down the mercury data (true values)
	ofp = open(caseID + '_' + 'hg_data.pkl','wb')
	pickle.dump(HG_reg,ofp)
	ofp.close()
	# save down loadings
	ofp = open(caseID + '_' + 'LOADINGS.pkl','wb')
	loads = np.array(HGresults[nam1.index('loadings')])
	pickle.dump(loads,ofp)
	ofp.close()
	# save R2 data on all components (NB--> default includes an intercept!!!!!)
	ofp = open(caseID + '_' + 'R2D2.pkl','wb')
	r2d2 = r.R2(HGresults)
	namr2 = list(r2d2.names)
	r2out = np.array(r2d2[namr2.index('val')]).squeeze()
	pickle.dump(r2out,ofp)
	ofp.close()
	# save the coefficients down.  (NB--> it's a matrix that is vars x ncomp)
	ofp = open(caseID + '_' + 'COEFFICIENTS.pkl','wb')
	allCoef = np.array(HGresults[nam1.index('coefficients')]).squeeze()
	pickle.dump(allCoef,ofp)
	ofp.close()
	# save down the X matrix
	ofp = open(caseID + '_' + 'Xmat.pkl','wb')
	pickle.dump(RHS,ofp)
	ofp.close()
	
		
	# make pngs of the predictions
	rootnm = 'HGreg_preds'
	if (cstudy < 9999):
		rootnm += ('_' + currstud)
	for i in xrange(nvars):
		preds_fn = newpath + rootnm + '_nc' + str(i+1) + '.png'
		r.png(file=preds_fn,height=750,width=750,)
		r.predplot(HGresults,asp=1,ncomp=i+1,line="True")
		r('''dev.off()''')
	
	predpng = 'HGreg_pred_scaled.png'
	predpng = newpath + predpng
	XYlims = (-3.5,1.5)
	#HGpred = np.array(list(r.fitted(HGresults))[len(HG_reg)*3:])
	HGpred = cv_data[:,opt_comp]
	mypredplot(5150,HG_reg,HGpred,'Predicted vs. Modeled: ' + str(opt_comp) + ' components',predpng,XYlims)
	
	
	
	
	# make png of RMSEP
	RMSEP_FN = 'HGreg_RMSEP'
	if (cstudy < 9999):
		RMSEP_FN += ('_' + currstud)
	RMSEP_FN += '.png'
	RMSEP_FN = newpath  + RMSEP_FN
	r.png(file=RMSEP_FN,height=750,width=750)
	r.plot(r.RMSEP(HGresults),legendpos="topright")
	r('''dev.off()''')
	
	
	a=r.seq(1,nvars)
	
	# look at scores
	#
	SCORE_FN = 'HGreg_SCORES'
	if (cstudy < 9999):
		SCORE_FN += ('_' + currstud)
	SCORE_FN += '.png'
	SCORE_FN = newpath  + SCORE_FN
	r.png(file=SCORE_FN,height=750,width=750)
	r.plot(HGresults,plottype="scores",comps=a)
	r('''dev.off()''')
	
	
	# look at loadings
	#
	LOAD_FN = 'HGreg_LOADINGS'
	if (cstudy < 9999):
		LOAD_FN += ('_' + currstud)
	LOAD_FN += '.png'
	LOAD_FN = newpath  + LOAD_FN
	r.png(file=LOAD_FN,height=750,width=750)
	r.plot(HGresults,plottype="loadings",comps=a,legend="bottomleft")
	r('''dev.off()''')
	

	
	
	#print out coefficients and r^2 values
	curdir = os.path.split(os.getcwd())[1]
	sumfile = curdir+ '_' + boundname + '_' + caseID + '_' + '_summary.dat'
	if (cstudy < 9999):
		sumfile = 'summary_' + currstud + '.dat'
	sumfile = newpath  + sumfile
	ofp = open(sumfile, 'w')
	ofp.write('summary and R^2 for ' + sumfile + '\n')
	ofp.write('\n'*4)
	ofp.write('*'*12)
	ofp.write('\nCOEFFICIENTS\n')
	ofp.write('%17s' %('VarName'))
	for i in np.arange(len(xname)):
	    ofp.write('%17s' %(str(i+1) + 'comps'))
	ofp.write('\n')
	for i in np.arange(len(xname)):
	    ofp.write('%17s' %(xname[i]))
	    for j in np.arange(len(allCoef[i,:])):
		ofp.write('%17f' %(allCoef[i,j]))
	    ofp.write('\n')
	ofp.write('*'*12)
	ofp.write('\nR^2\n')
	for i,v in enumerate(r2out):
	    ofp.write('%-10s-->' %(str(i) + ' comps'))
	    ofp.write('%10f\n' %(v))
	ofp.close()
	
	
	
	#######################
	#   MAKE PREDICTIONS  #
	#######################
	
	if PredFlag:
		# first load up the prediction data
		indata = np.loadtxt(pred_infile, skiprows=1,delimiter='~')
		pHpr    = np.atleast_2d(indata[:,1])
		so4pr   = np.atleast_2d(indata[:,2])
		TOCpr   = np.atleast_2d(indata[:,3])
		WETpr   = np.atleast_2d(indata[:,4])
	
		if so4logFlag:
			so4pr = np.log10(so4pr)
			so4pr[np.isinf(so4pr)] = -1
		if toclogFlag:
			TOCpr = np.log10(TOCpr)
			TOCpr[np.isinf(TOCpr)] = -2
		if wetlandlogFlag:
			WETpr = np.log10(WETpr)
			WETpr[np.isinf(WETpr)] = -4
		
		if wetlflag:
			RHSpr = np.hstack((so4pr.T,TOCpr.T,pHpr.T,WETpr.T))
		else:
			RHSpr = np.hstack((so4pr.T,TOCpr.T,pHpr.T))
		junkus = np.zeros_like(so4pr.T)       
		RHSpr    = rob.conversion.py2ri(RHSpr)
		junkus   = rob.conversion.py2ri(junkus)
		
		PredData = rlc.TaggedList([junkus,RHSpr],tags=('junkus','rhs'))
		PredData = rob.DataFrame(PredData)
		rob.globalenv["PredData"] = PredData
		rpred = np.array(r.predict(HGresults, ncomp=opt_comp, newdata=PredData)).squeeze()
		ofp = open(newpath + curdir + '_' + boundname  + '_' + caseID + '_' + 'predvals.dat', 'w')
		ofp.write('stationID~predicted_log10_methyl_HG\n')
		for i,line in enumerate(rpred):
			ofp.write(str(indata[i,0]) + "~" + str(line) + '\n')
		ofp.close()