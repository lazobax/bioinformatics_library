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

  return PAM250[xi][yi]

def visualizeMatrix(matrix):

  rowLenght = len(matrix[0])
  columnLength = len(matrix)

  for col in range(columnLength):
    print(matrix[col])

def initializeMatrix(Seq1, Seq2, PENALTY):

  matrix = [[None for i in range(len(Seq1))] for j in range(len(Seq2))]

  matrix[0][0] = 0

  for i in range(1, len(matrix)):

    row = matrix[i]

    row[0] = PENALTY*i

  for i in range(1, len(matrix[0])):

    row1 = matrix[0]

    row1[i] = PENALTY*i

  return matrix

def GlobalAlignment(Seq1, Seq2, PENALTY):

  matrix = initializeMatrix(Seq1, Seq2, PENALTY)
  tracker = [[None for i in range(len(Seq1))] for j in range(len(Seq2))]

  optimalCoordinates = [0,0]
  optimal = 0

  for j in range(1, len(matrix)):
    for i in range(1, len(matrix[0])):
      diagonal = matrix[j-1][i-1] + PAM250(Seq1[i], Seq2[j])
      right = matrix[j][i-1] + PENALTY
      down = matrix[j-1][i] + PENALTY

      choice = max(diagonal, right, down)

      if choice == diagonal:
        tracker[j][i] = 'Diagonal'
      elif choice == right:
        tracker[j][i] = 'Right'
      else:
        tracker[j][i] = 'Down'

      if(choice > optimal):
        optimal = choice
        optimalCoordinates = [j, i]

      matrix[j][i] = choice


  aligned1 = ""
  aligned2 = ""
  cD = tracker[optimalCoordinates[0]][optimalCoordinates[1]]
  cC = optimalCoordinates

  while cD != None:
    if cD == 'Diagonal':
      aligned1 += Seq1[cC[1]]
      aligned2 += Seq2[cC[0]]
      cC = [cC[0]-1, cC[1]-1]
    elif cD == 'Right':
      aligned1 += Seq1[cC[1]]
      aligned2 += "-"
      cC = [cC[0], cC[1]-1]
    else:
      aligned1 += "-"
      aligned2 += Seq2[cC[0]]
      cC = [cC[0]-1, cC[1]]
    
    cD = tracker[cC[0]][cC[1]]

  print(aligned1[::-1])
  print(aligned2[::-1])

Seq1 = "-ATWES" 
Seq2 = "-TCAET"
PENALTY = -2

GlobalAlignment(Seq1, Seq2, PENALTY)








