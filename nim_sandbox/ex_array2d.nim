type
  Array2D[T] = object
    xmin, xmax: int
    ymin, ymax: int
    ny: int           # size of 2nd dimension
    data: seq[T]      # contiguous storage

proc initArray2D[T](xmin, xmax, ymin, ymax: int): Array2D[T] =
  result.xmin = xmin
  result.xmax = xmax
  result.ymin = ymin
  result.ymax = ymax
  result.ny = ymax - ymin + 1
  let nx = xmax - xmin + 1
  result.data = newSeq[T](nx * result.ny)

proc idx[T](A: Array2D[T], x, y: int): int =
  assert x >= A.xmin and x <= A.xmax
  assert y >= A.ymin and y <= A.ymax
  (x - A.xmin) * A.ny + (y - A.ymin)

proc `[]`[T](A: Array2D[T], x, y: int): T =
  A.data[A.idx(x,y)]

proc `[]=`[T](A: var Array2D[T], x, y: int, value: T) =
  A.data[A.idx(x,y)] = value


var A = initArray2D[float](-5, 5, 3, 7)

A[-2,4] = 3.14
echo A[-2,4]

# Fill entire grid
for i in -5..5:
  for j in 3..7:
    A[i,j] = float(i*j)

echo A[0,5]


# Vector of offset arrays
#var grids: seq[Array2D[float]]
#var grids = newSeq

var grids: seq[Array2D[float]] = @[]
grids.setLen(3)

for k in 0..<3:
  grids[k] = initArray2D[float](-5,5,-5,5)

grids[0][-2,3] = 1.0
echo grids[0][-2,3]
echo grids[1][-2,3]




