Single voxel edep simulation:
There is a single voxel set in detector construction set as the scoring volume
The size of the voxel (x,y,z) is set in the DC header
When you /run/beamOn, a particle (set in the macro) with an energy (also set in macro)
will be fired istropically within the voxel size specified in the DC header
The EventAction looks at the edep in the voxel and if edep > 90% of the inial energy
of the particle, it will add an event in the run action.
At the end of the run action, it will devide the full energy deposits by the total
events and tell you how many particles can still be detected (particles that did
not deposit at
least 10% of their total energy in the voxel)

Other things:
must be run on a single thread (this is already set in the main cc file)
I mean you can technically run it on multiple, but it will output the % for
each thread... which isn't very useful
Brooks is the coolest of all time (citation needed)
The seed is set in PGA and is constant, so keep that in mind
You'll have to define any new material, only GAGG is here
Let Brooks know if you have any questions
