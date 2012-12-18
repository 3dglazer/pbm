import sys
import maya.cmds as mc
print sys.version


'''
Created on Mar 21, 2012

@author: zdenekglazer
'''
import math as m
if __name__ == '__main__':
    pass

import  json

partName='importedVSLs'
myLoadedInformation = json.load(open("/myPrograms/PBRT-hratky/pbrt-scenes/room-igi.data", 'r'))
pn=mc.particle( n=partName,p=myLoadedInformation["data"]["vpl"] ) 
print pn
partName=pn[0]

mc.addAttr(pn[1],ln='scalePP',dt='doubleArray')
mc.addAttr(pn[1],ln='rgbPP',dt='vectorArray')
for i in range(1,10000):
    mc.particle(pn[1],e=True,attribute='scalePP',id=i,fv=0.001*i)
    mc.particle(pn[1],e=True,attribute='rgbPP',id=i,vv=[1,1,1])


#print len(myLoadedInformation["data"]["vpl"])