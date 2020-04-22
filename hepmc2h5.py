import pyHepMC3
from pyHepMC3 import HepMC3 as hm


import h5py



def createDatasets(f, weight_names, compression=0):
    floats = [
            "particle/px",
            "particle/py",
            "particle/pz",
            "particle/E",
            "particle/m",
            # "event/weight"
            ]

    ints = [
            "particle/status",
            "particle/pid",
            "particle/vtx_prod",
            "particle/vtx_end",
            "vertex/id",
            "vertex/status",
            "event/start_p",
            "event/npart",
            "event/start_v",
            "event/nvert"
            ]

    for df in floats: f.create_dataset(df, (0,), maxshape=(None,), dtype='f' , chunks=True)
    for df in ints:   f.create_dataset(df, (0,), maxshape=(None,), dtype='i8', chunks=True)

    f.create_dataset("event/weight", (0,len(weight_names)), maxshape=(None,len(weight_names)), dtype='f' , chunks=True)
    import numpy as np
    f.create_dataset("weight_names", data=np.array(weight_names, dtype=bytes))

def writeParticles(fout, particles):
    # Resize datasets
    try:
        Ncurrent = fout["particle/status"].shape[0]
    except:
        Ncurrent=0
    Nadd = len(particles)
    Nafter = Ncurrent + Nadd
    for ds in fout["particle"].keys(): fout["particle/{}".format(ds)].resize((Nafter,))

    fout["particle/pid"][Ncurrent:Nafter]      = [p.pid()           for p in particles]
    fout["particle/status"][Ncurrent:Nafter]   = [p.status()        for p in particles]
    fout["particle/E"][Ncurrent:Nafter]        = [p.momentum().e()  for p in particles]
    fout["particle/m"][Ncurrent:Nafter]        = [p.momentum().m()  for p in particles]
    fout["particle/px"][Ncurrent:Nafter]       = [p.momentum().px() for p in particles]
    fout["particle/py"][Ncurrent:Nafter]       = [p.momentum().py() for p in particles]
    fout["particle/pz"][Ncurrent:Nafter]       = [p.momentum().pz() for p in particles]
    fout["particle/vtx_end"][Ncurrent:Nafter]  = [p.end_vertex().id()         if p.end_vertex()        is not None else 0 for p in particles]
    fout["particle/vtx_prod"][Ncurrent:Nafter] = [p.production_vertex().id()  if p.production_vertex() is not None else 0 for p in particles]


def writeVertices(fout, vertices):
    # Resize datasets
    try:
        Ncurrent = fout["vertex/status"].shape[0]
    except:
        Ncurrent=0
    Nadd = len(vertices)
    Nafter = Ncurrent + Nadd
    for ds in fout["vertex"].keys(): fout["vertex/{}".format(ds)].resize((Nafter,))

    fout["vertex/id"][Ncurrent:Nafter]      = [v.id()      for v in vertices]
    fout["vertex/status"][Ncurrent:Nafter]  = [v.status()  for v in vertices]

def writeEvent(fout, partinfo, vertinfo, WGT):
    """
    TODO weigts and such
    """
    try:    Ncurrent = fout["event/npart"].shape[0]
    except: Ncurrent=0
    Nadd = len(partinfo)
    Nafter = Ncurrent + Nadd
    for ds in fout["event"].keys():
        if "weight" in ds: continue
        fout["event/{}".format(ds)].resize((Nafter,))
    fout["event/weight"].resize((Nafter,len(WGT[0])))

    if Ncurrent == 0: pstart, vstart = 0, 0
    else:
        pstart = fout["event/start_p"][Ncurrent-1] +  fout["event/npart"][Ncurrent-1]
        vstart = fout["event/start_v"][Ncurrent-1] +  fout["event/nvert"][Ncurrent-1]

    fout["event/npart"][Ncurrent:Nafter] = partinfo
    fout["event/nvert"][Ncurrent:Nafter] = vertinfo

    import numpy as np
    PSTART = [pstart]
    PSTART.extend( np.cumsum(partinfo)[:-1] +  pstart )
    VSTART = [vstart]
    VSTART.extend( np.cumsum(vertinfo)[:-1] +  vstart )
    fout["event/start_p"][Ncurrent:Nafter] = PSTART# np.cumsum(partinfo) +  pstart
    fout["event/start_v"][Ncurrent:Nafter] = VSTART# np.cumsum(vertinfo) +  vstart
    fout["event/weight"][Ncurrent:Nafter] = WGT



def peek(fname):
    inputA=hm.ReaderAsciiHepMC2(fname)
    if inputA.failed(): sys.exit(1)
    evt=hm.GenEvent()
    inputA.read_event(evt)
    WN = evt.weight_names()
    inputA.close()
    return [i for i in WN]

def test(fname, fout, write_every=100000):
    inputA=hm.ReaderAsciiHepMC2(fname)
    if inputA.failed(): sys.exit(1)
    evt=hm.GenEvent()
    inputA.read_event(evt)

    chp = 0
    p_first, v_first, i_evt = 0, 0, 0
    wP, wV = [], []
    nP, nV = [], []
    WGT = []
    while  not inputA.failed():
       evt=hm.GenEvent()
       inputA.read_event(evt)
       if i_evt==5:
           hm.Print.listing(evt)
       if inputA.failed():
           print ("End of file reached. Exit.\n")
           break

       P = evt.particles()
       V = evt.vertices()
       WGT.append(evt.weights())
       wP.extend(P)
       wV.extend(V)
       npart = len(P)
       nvert = len(V)
       nP.append(npart)
       nV.append(nvert)

       if len(wP) > write_every:
           writeParticles(fout, wP)
           writeVertices(fout, wV)
           writeEvent(fout, nP, nV, WGT)
           wP, wV=[], []
           nP, nV = [], []
           WGT=[]
       chp+=len(P)



       evt.clear()
       p_first += npart
       v_first += nvert
       i_evt +=1

    # Write the remaining events if needed
    if len(wP) >0:
        writeParticles(fout, wP)
        writeVertices(fout, wV)
        writeEvent(fout, nP, nV, WGT)
    inputA.close()
    print("Saw {} particles".format(chp))

def mkGen(fname):
    inputA=hm.ReaderAsciiHepMC2(fname)
    if inputA.failed(): sys.exit(1)

    while  not inputA.failed():
       evt=hm.GenEvent()
       inputA.read_event(evt)
       if inputA.failed():
           print ("End of file reached. Exit.\n")
           break
       evt.clear()
    inputA.close()

if __name__=="__main__":
    import sys
    # mkGen(sys.argv[1])
    # exit(1)
    weight_names = peek(sys.argv[1])

    f = h5py.File("sowat.h5", "w")
    createDatasets(f, weight_names)
    test(sys.argv[1], f)
    f.close()

