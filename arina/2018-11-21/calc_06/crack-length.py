import medcoupling as mc

medFileName = 'CT_results.rmed'
fnout       = 'crack-length6.dat'
nume        = (1,)
displName   = 'RESU____DEPL'
fieldName   = 'RESU____VARI_ELGA'
numCmp      = 1                         # porosity
seuil       = 0.2
xtip        = 15.2
ytip        = 0.0
small       = 1.e-10


meshName   = mc.GetMeshNamesOnField(medFileName,fieldName)[0]
timeStamps = mc.GetAllFieldIterations(medFileName,fieldName)


for ts in timeStamps:
    first = (ts[0] == timeStamps[0][0])
    
    #if ts[0] not in nume: continue
    #ts = timeStamps[900]
    print ts
    
    
    # Champ d'endommagement
    vari = mc.ReadField(mc.ON_GAUSS_PT,medFileName,meshName,0,fieldName,ts[0],ts[1]).keepSelectedComponents([numCmp])
    poromax = vari.getArray().getMaxValue()[0]


    # Id du noeud tip
    if first:
        xy = vari.getMesh().getCoords()
        idx = xy[:,0].findIdsInRange(xtip-small,xtip+small)
        idy = xy[:,1].findIdsInRange(ytip-small,ytip+small)
        idt = idx.buildIntersection(idy)
        assert len(idt)==1
    
    
    if poromax < seuil: 
        da = 0.0
        
    else:
    
        # Champ de deplacement
        displField = mc.ReadField(mc.ON_NODES,medFileName,meshName,0,displName,ts[0],ts[1])
        displField.setNature(mc.IntensiveMaximum)


        # Projection du deplacement sur le maillage du champ d'endommagement
        srcMesh = displField.getMesh()
        tarMesh = vari.getMesh()
        srcMesh.convertQuadraticCellsToLinear()
        srcMesh.simplexize(0)
        displ = displField.getValueOnMulti(tarMesh.getCoords())
        u = displ.keepSelectedComponents((0,1))

        #remap = mc.MEDCouplingRemapper()
        #remap.prepare(srcMesh,tarMesh,"P1P1")
        #displ = remap.transferField(displField,1e300)
        #u = displ.getArray().keepSelectedComponents((0,1))


        # Actualisation du maillage par le champ de deplacement
        mesh = tarMesh
        xy = mesh.getCoords()
        u.setInfoOnComponents(xy.getInfoOnComponents())
        xy += u


        # Position de la pointe
        a0 = xy[idt[0]][0,0]


        # Champ constant pas cellule autour des points de Gauss
        poro  = vari.voronoize(1e-12)
        gmsh  = poro.getMesh()
        ctr   = gmsh.computeIsoBarycenterOfNodesPerCell()    


        # Cellules correspondant a des points casses
        ids = poro.findIdsInRange(seuil, 1.0)
        if len(ids) == 0:
            a = a0
        else:
            zone = ctr[ids]
            x = zone[:,0]
            a = x.getMaxValue()[0]
        da = a-a0
        
        
    if first:
        mode = 'w'
    else:
        mode = 'a'
    with file(fnout,mode) as fl:
        fl.write("%22.15f  %22.15f \n" % (ts[2],da))


    
    
    #vc.setNature(mc.IntensiveMaximum)
    
    ## Maillage P1 pour la destination
    #mesh.convertQuadraticCellsToLinear()
    
    ## Projecteur
    #if not projec:
        #projec = True
        #print 'Construction du projecteur...',
        #remap = rmp.MEDCouplingRemapper()
        #remap.prepare(vc.getMesh(),mesh,"P0P1")
        #print 'finie'
        
    #vn = remap.transferField(vc,0.0)
    #vn.setName('damage')
    #vn.getArray().setInfoOnComponent(0,'damage')

    #fileName = outFileBase + '_nd_' + repr(ts[0]) + '.med'
    #mc.WriteField(fileName,vn,True)
    #print ts[0],ts[2]


