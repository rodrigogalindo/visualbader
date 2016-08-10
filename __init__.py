## Extension to read a file output from AIMALL
## based on the code ReadXYZ
def readAIM(fileName):
        import re
	from OpenSave import osOpen
	from chimera import UserError, Molecule, Element, openModels, replyobj
	from chimera import Coord, connectMolecule, MaterialColor
        anums = {}

        ## Open the file
	getatoms = osOpen(fileName)

        ## Search for the string:
        ## 'Number of NACPs'
        ## which gives us the number of atoms
        ## and store it in numAtoms
        for data in getatoms:
                data = data.strip()
                if data.startswith("Number of NACPs"):
                        a,b,c,d,e = data.split()
                        numAtoms = int(e)
        getatoms.close()

        ###################################################################
	f = osOpen(fileName)
        readAtoms = 0
        state = "init"

        ## Skip the 26 header lines from the mpgviz file
        #f = f.readlines()[4:]

        for lines in f:
                lines = lines.strip()
                if lines.startswith("Atom      Charge                X                  Y                  Z"):
                        f.next()
                        for line in f:
                                if readAtoms == numAtoms:
                                        break
                                line = line.strip()       
                                if state == "init":
                                        state = "post"
                                        m = Molecule()
                                        r = m.newResidue("UNK", " ", 1, " ")
                                        state = "atoms"
                                        serial = 1
                                if not line: continue

                                elem, charge, x, y, z = line.split()
                                x, y, z = [float(c) for c in [x,y,z]]
                                ## Remove the number from the atomic symbol
                                elem = filter(lambda x: x.isalpha(), elem)
                                element = Element(elem)
                                anum = anums.get(element.name, 0) + 1
                                anums[element.name] = anum
                                a = m.newAtom("%s%d" % (element.name, anum), element)
                                r.addAtom(a)
                                a.setCoord(Coord((x*0.52918), (y*0.52918), (z*0.52918)))
                                a.serialNumber = serial
                                serial += 1
                                readAtoms += 1
                        connectMolecule(m)

        
        ###################################################################
        ## Critical points
        idNACP = 1
        mNACP = Molecule()
        rNACP = mNACP.newResidue("Nuclear atractors critical points", " ", 1, " ")

        idBCP = 1
        mBCP = Molecule()
        rBCP = mBCP.newResidue("Bond critical points", " ", 1, " ")

        idBPP = 1
        mBPP = Molecule()
        rBPP = mBPP.newResidue("Bond path points", " ", 1, " ")

        cpnumber = 1
        cpread = 0
	getcriticalpoints = osOpen(fileName)
        for text in getcriticalpoints:
                text = text.strip()
                if cpread == 0:
                        if text.startswith("CP#"):
                                readtext = text
                                cpread = 1
                elif cpread == 1:
                        if text.startswith("Type = (3,-3)"):
                                a,b,signature,type,atom1 = text.split()
                                cpnum,id,a,b,x,y,z = readtext.split()
                                x, y, z = [float(c) for c in [x,y,z]]
                                
                                ## Remove the number from the atomic symbol
                                elem = filter(lambda x: x.isalpha(), atom1)
                                element = Element(elem)
                                anum = anums.get(element.name, 0) + 1
                                anums[element.name] = anum
                                #atomNACP = m.newAtom("%s%d" % (element.name, anum), element)
                                atomNACP = mNACP.newAtom("%s%s" % ("NACP# ", atom1), element)
                                atomNACP.serialNumber = serial
                                #atomNACP = mNACP.newAtom(atom1, Element("NACP"))
                                atomNACP.drawMode = 1
                                atomNACP.radius = 0.2
                                rNACP.addAtom(atomNACP)
                                atomNACP.setCoord(Coord((x*0.52918), (y*0.52918), (z*0.52918)))
				cpread = 0
				serial += 1
                        elif text.startswith("Type = (3,+3)"):
                                print 'CCP'
                                cpread = 0
                        elif text.startswith("Type = (3,+1)"):
                                print 'RCP'
                                cpread = 0
                        elif text.startswith("Type = (3,-1)"):
                                a,b,signature,type,atom1,atom2 = text.split()
                                cpnum,id,a,b,x,y,z = readtext.split()
                                x, y, z = [float(c) for c in [x,y,z]]                                

                                ## Remove the number from the atomic symbol
                                elem = filter(lambda x: x.isalpha(), atom1)
                                element = Element(elem)
                                anum = anums.get(element.name, 0) + 1
                                anums[element.name] = anum
				atomBCP = mBCP.newAtom("%s%s" % ("BCP# ", id), element)
                                atomBCP.serialNumber = serial
                                #atomBCP = mBCP.newAtom(atom1, Element("He"))
                                atomBCP.color = MaterialColor(1,0,0,1)
				atomBCP.drawMode = 1
                                atomBCP.radius = 0.07
                                rBCP.addAtom(atomBCP)
                                atomBCP.setCoord(Coord((x*0.52918), (y*0.52918), (z*0.52918)))
                                cpread = 0
                                serial += 1

        ###################################################################
        ## Bond Paths

        getbondpaths = osOpen(fileName)
        statusSamples = 0
        countSamples = 0
	totalSamples = 0
        skip = 0
	for samples in getbondpaths:
                samples = samples.strip()
                if re.match("(.*)sample points along path from BCP to atom(.*)", samples):
                        numSamples,a,b,c,d,e,tipo,g,h,atom = samples.split()
                        numSamples = int(numSamples)
                        statusSamples = 1
			totalSamples = int(totalSamples)

                if statusSamples == 1:
			totalSamples = totalSamples + statusSamples
                        if countSamples == numSamples:
                                countSamples = 0
                                statusSamples = 0
                        else:
                                try:
                                        x,y,z,rho = samples.split()
                                        x, y, z = [float(c) for c in [x,y,z]]
					
					
                                        ## Remove the number from the atomic symbol
					elem = filter(lambda x: x.isalpha(), atom1)
					element = Element(elem)
					anum = anums.get(element.name, 0) + 1
					anums[element.name] = anum
					atomBPP = mBPP.newAtom("%s%s" % ("BPP", serial), element)
					atomBPP.serialNumber = serial
                                        #atomBPP = mBPP.newAtom(atom1, Element("kr"))
                                        #atomBPP = mBPP.newAtom("%s%s" % ("BPP", serial), element)
					atomBPP.color = MaterialColor(0,206,209,1)
					atomBPP.drawMode = 1
					atomBPP.radius = 0.01
					rBPP.addAtom(atomBPP)
					atomBPP.setCoord(Coord((x*0.52918), (y*0.52918), (z*0.52918)))
					#cpread = 0
					serial += 1
					#cpnumber += 1
					countSamples += 1
				except ValueError:
                                        countSamples += 1
                                        continue


        mNACP.isRealMolecule = False
#        mNACP.noprefs = True
        mBCP.isRealMolecule = False
        mBCP.noprefs = True
        mBPP.isRealMolecule = False
        mBPP.noprefs = True
      

        ## Return molecules
        return [m, mNACP, mBCP, mBPP]
