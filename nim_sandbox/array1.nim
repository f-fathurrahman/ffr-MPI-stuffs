# 2D array with arbitrary starting indices
var A: array[-2..2, array[3..5, float]]

# Fill the array
for i in -2..2:
  for j in 3..5:
    A[i][j] = float(i * j)

# Access example
echo A[-2][3]
echo A[0][4]


