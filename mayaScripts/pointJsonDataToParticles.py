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


myLoadedInformation = json.load(open("/myPrograms/PBRT-hratky/pbrt-scenes/room-igi.data", 'r'))
mc.particle( p=myLoadedInformation["data"]["vpl"] ) 
print len(myLoadedInformation["data"]["vpl"])