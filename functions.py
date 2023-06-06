

#reads FASTA file and returns a list of sequences
def readFastaFile(fileName):

  fileObj = open(fileName, 'rU')
  sequences = []
  seqFragments = []

  for line in fileObj:
    #if line starts with >, it is a header line
    if line.startswith('>'):
      #if there is a sequence in seqFragments, join the fragments together and add to sequences and reset the fragments
      if len(seqFragments) > 0:
        sequences.append(''.join(seqFragments))
        seqFragments = []
    #if line does not start with >, it is a sequence fragment
    else:
      seqFragments.append(line.rstrip()) #rstrip removes trailing whitespace

  #add the last sequence to sequences
  if len(seqFragments) > 0:
    sequences.append(''.join(seqFragments))

  fileObj.close()

  return sequences

#reads a PDB file and calculates the average x, y, z coordinates of the atoms in the protein
def calculateCentroid(pdbFile):
  fileObj = open(pdbFile, 'rU')

  nAtoms = xsum = ysum = zsum = 0  #initialize variables to 0

  for line in fileObj:
    if line[:6] == 'ATOM  ':#ATOM is the first 6 characters of the line
      nAtoms += 1
      xsum += float(line[30:38]) #x coordinate is characters 31-38
      ysum += float(line[38:46]) #y coordinate is characters 39-46
      zsum += float(line[46:54]) #z coordinate is characters 47-54
    fileObj.close()

    if nAtoms == 0:
      xavg = yavg = zavg = 0
    else:
      xavg = xsum / nAtoms
      yavg = ysum / nAtoms
      zavg = zsum / nAtoms
      
  return xavg, yavg, zavg

def PAM250(x, y):
  
  AAs = ["C", "S", "T", "P", "A", "G", "N", 
          "D", "E", "Q", "H", "R", "K", "M",
          "I", "L", "V", "F", "Y", "W"]

  AminoHash = {AAs[i] : i for i in range(len(AAs))}

  PAM250 = [[12],
            [0,2],
            [-2,1,3],
            [-3,1,0,6],
            [-2,1,1,1,2],
            [-3,1,0,-1,1,5],
            [-4,1,0,-1,0,0,2],
            [-5,0,0,-1,0,1,2,4],
            [-5,0,0,-1,0,0,1,3,4],
            [-5,-1,-1,0,0,-1,1,2,2,4],
            [-3,-1,-1,0,-1,-2,2,1,1,3,6],
            [-4,0,-1,0,-2,-3,0,-1,-1,1,2,6],
            [-5,0,0,-1,-1,-2,1,0,0,1,0,3,5],
            [-5,-2,-1,-2,-1,-3,-2,-3,2,-1,-2,0,0,6],
            [-2,-1,0,-2,-1,-3,-2,-2,-2,-2,-2,-2,-2,2,5],
            [-6,-3,-2,-3,-2,-4,-3,-4,-3,-2,-2,-3,-3,4,2,6],
            [-2,-1,0,-1,0,-1,-2,-2,-2,-2,-2,-2,-2,2,4,2,4],
            [-4,-3,-3,-5,-4,-5,-4,-6,-5,-5,-2,-4,-5,0,1,2,-1,9],
            [0,-3,-3,-5,3,-5,-2,-4,-4,-4,0,-4,-4,-2,-1,-1,-2,7,10,],
            [-8,-2,-5,-6,-6,-7,-4,-7,-7,-5,-3,2,-3,-4,-5,-2,-6,0,0,17]]

  xi = AminoHash[x]
  yi = AminoHash[y]

  #ensure that xi is greater than yi
  if (yi > xi):
    xi, yi = yi, xi

  return PAM250[yi][xi]

