scf_inp_template = '''###############################################################################
#  SPR-KKR input file    *SCF.inp 
###############################################################################

CONTROL  DATASET     = {system} 
         ADSI        = SCF 
         POTFIL      = {system}.pot 
         PRINT = 0    

MODE     SP-SREL 

TAU      BZINT= POINTS  NKTAB= 4000 

ENERGY   GRID={5}  NE={30} 
         ImE=0.0 Ry   EMIN=-0.2 Ry

SCF      NITER=200 MIX=0.05 VXC=PBE
         TOL=0.00001  MIXOP=0.20  ISTBRY=1 
         QIONSCL=0.80 
         NOSSITER
         MSPIN={magmoms} 
'''

replaces_scf = {'{system}': '',
                '{magmoms}': ''}

jxc_inp_template = '''###############################################################################
#  SPR-KKR input file    *JXC.inp 
#  created by xband on Sat 14 Jun 22:38:43 BST 2025
###############################################################################
 
CONTROL  DATASET     = {system} 
         ADSI        = JXC 
         POTFIL      = {system}.pot 
         PRINT = 0    
 
MODE     SP-SREL 
 
TAU      BZINT= POINTS  NKTAB= 2000 
 
ENERGY   GRID={5}  NE={30} 
         EMIN=-0.2 Ry
 
TASK     JXC   CLURAD={cluster_radius}
'''

replaces_jxc = {'{system}': '',
                '{cluster_radius}': ''}