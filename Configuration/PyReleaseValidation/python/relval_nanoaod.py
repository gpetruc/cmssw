# import the definition of the steps and input files:
from  Configuration.PyReleaseValidation.relval_steps import *

# here only define the workflows as a combination of the steps defined above:
workflows = Matrix()

# each workflow defines a name and a list of steps to be done. 
# if no explicit name/label given for the workflow (first arg),
# the name of step1 will be used

workflows[70000801] = ['', ['TT_2016_80X_NANO_INPUT','NANOAODMC2016_80X']]
workflows[70000921] = ['', ['TT_2017_92X_NANO_INPUT','NANOAODMC2017_92X']]
workflows[70000941] = ['', ['TT_2017_94X_NANO_INPUT','NANOAODMC2017'    ]]

workflows[70100801] = ['', ['RunJetHT2016H_NANO_INPUT','NANOAOD2016_80X']]
workflows[70100802] = ['', ['RunMET2016H_NANO_INPUT','NANOAOD2016_80X']]
workflows[70100803] = ['', ['RunSingleEl2016H_NANO_INPUT','NANOAOD2016_80X']]
workflows[70100804] = ['', ['RunSingleMu2016H_NANO_INPUT','NANOAOD2016_80X']]
workflows[70100921] = ['', ['RunJetHT2017C_R299649_NANO_INPUT','NANOAOD2017_92X']]
workflows[70100922] = ['', ['RunMET2017C_R299649_NANO_INPUT','NANOAOD2017_92X']]
workflows[70100923] = ['', ['RunSingleEl2017C_R299649_NANO_INPUT','NANOAOD2017_92X']]
workflows[70100924] = ['', ['RunSingleMu2017C_R299649_NANO_INPUT','NANOAOD2017_92X']]

