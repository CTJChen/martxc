#### 2019/07/11 
Version 0.1
Passed initial tests on an OSX laptop and a Fedora 29 linux desktop,
with both python 2.7 and 3.6.8

2019/07/12 - V 0.1.1 
martxcexpmap.py can now generate a exposure time weighed off-axis angle map. 
martmkarf.py now recognize the off-axis map by martxcexpmap, and the off-axis and aperture corrections can be based on the time-averaged off-axis angle map. This primarily for the all-sky survey observations.