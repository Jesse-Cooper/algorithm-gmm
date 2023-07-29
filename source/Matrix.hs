-- Purpose  Provides a vector and matrix data structure with functions to
--          manipulate them.


module Matrix (Vector, Matrix,
               validMatrix,
               nRows, nCols,
               mapMatrix,
               mapMatrixWith,
               mapDiagonals,
               foldRows, foldCols,
               zipMatrixWith, zipRowsWith, zipColsWith,
               det2x2, inv2x2, cholesky2x2,
               outerProd,
               quadForm) where


-- Both `Vector` and `Matrix` are just lists.

-- [[a, b],   * [a, b] and [c, d] are separate rows
--  [c, d]]   * a/c and b/d are in the same column

-- A valid matrix is defined as having the same number of columns in each row.
-- `validMatrix` provides a way to check if a matrix is valid.
type Vector a = [a]
type Matrix a = Vector (Vector a)


-- | Determines if a matrix is valid.
--
--   Checks each row has the same length (number of columns) as the first row.
validMatrix :: Matrix a -> Bool
validMatrix matrix = all ((nCols matrix ==) . length) matrix


-- | Finds the number of rows in a matrix.
--
--   Matrices are organised by rows, so the length of the outer list is the
--   number of rows.
nRows :: Matrix a -> Int
nRows matrix = length matrix


-- | Finds the number of columns in a matrix.
--
--   Matrices are organised by rows, all with the same length, so the length of
--   the first row gives the number of columns.
nCols :: Matrix a -> Int
nCols [] = 0
nCols (row:_) = length row


-- | Applies a function to a single element in a vector by its index.
--
--   Linearly traverses list until at the element (when `i` is 0), then applies
--   the function.
--   `i` represents the number of elements remaining to traverse.
mapVecAtIndex :: (a -> a) -> Vector a -> Int -> Vector a
mapVecAtIndex _ [] _ = error "out of index"
mapVecAtIndex f (x:xs) i
    | i <= 0    = f x : xs
    | otherwise = x : mapVecAtIndex f xs (i - 1)


-- | Applies a function to each element in a matrix.
--
--   A double map touches each element as a matrix is a list of lists.
mapMatrix :: (a -> b) -> Matrix a -> Matrix b
mapMatrix f matrix = map (map f) matrix


-- | Applies a function with a constant to each element in a matrix.
--
--   Just uses `mapMatrix` but currys the given function with the constant.
--
--   mapMatrixWith (+) matrix int
--   [[a, b], + 1 = [[a + 1, b + 1],
--    [c, d]]        [c + 1, d + 1]]
mapMatrixWith :: (a -> b -> c) -> Matrix a -> b -> Matrix c
mapMatrixWith f matrix x = mapMatrix (flip f x) matrix


-- | Applies a function to the diagonals of a matrix.
--
--   Diagonals are elements where the index of the row == index of the column.
--   Linearly traverses down each row and applies the function to the element of
--   the row's index in the outer matrix list.
mapDiagonals :: (a -> a) -> Matrix a -> Matrix a
mapDiagonals f matrix
    | nRows matrix /= nCols matrix
        = error "matrix must be square"
    | otherwise
        = zipWith (mapVecAtIndex f) matrix [0..]


-- | Folds each row of a matrix into a single vector.
--
--   Maps a fold over each row of the matrix.
--
--   foldRows (+) matrix
--   [[a, b], = [a + b, c + d]
--    [c, d]]
foldRows :: (a -> a -> a) -> Matrix a -> Vector a
foldRows f matrix = map (foldl1 f) matrix


-- | Folds each column of a matrix into a single vector.
--
--   Folds the entire matrix with a zip of the function where the zip matches
--   the columns to each row.
--
--   foldCols (+) matrix
--   [[a, b], = [a + c, b + d]
--    [c, d]]
foldCols :: (a -> a -> a) -> Matrix a -> Vector a
foldCols f matrix = foldl1 (zipWith f) matrix


-- | Element-wise zips two matrices of the same size with a function.
--
--   Zips both matrices by their rows then zips their columns within each row.
--
--   zipMatrixWith (+) matrix matrix
--   [[a, b], + [[1, 2], = [[a + 1, b + 2],
--    [c, d]]    [3, 4]]    [c + 3, d + 4]]
zipMatrixWith :: (a -> b -> c) -> Matrix a -> Matrix b -> Matrix c
zipMatrixWith f x y
    | nRows x /= nRows y ||  nCols x /= nCols y
        = error "matrices must be the same size"
    | otherwise
        = zipWith (zipWith f) x y


-- | Zips each row with the same vector with a function.
--
--   Maps each row with a zip of the vector. The vector must be the same length
--   as the number of columns in the matrix.
--
--   zipRowsWith (+) matrix vector
--   [[a, b], + [1, 2] = [[a + 1, b + 2],
--    [c, d]]             [c + 1, d + 2]]
zipRowsWith :: (a -> b -> c) -> Matrix a -> Vector b -> Matrix c
zipRowsWith f matrix vector
    | nCols matrix /= length vector
        = error "vector must be the same length as columns in matrix"
    | otherwise
        = map (\xs -> zipWith f xs vector) matrix


-- | Zips each column with the same vector with a function.
--
--   Matches the vector elements to the matrix rows by a zip, then maps those
--   elements across the rows. The vector must be the same length as the number
--   of rows in the matrix.
--
--   zipColsWith (+) matrix vector
--   [[a, b], + [1, 2] = [[a + 1, b + 1],
--    [c, d]]             [c + 2, d + 2]]
zipColsWith :: (a -> b -> c) -> Matrix a -> Vector b -> Matrix c
zipColsWith f matrix vector
    | nRows matrix /= length vector
        = error "vector must be the same length as rows in matrix"
    | otherwise
        = zipWith (\xs y -> map (flip f y) xs) matrix vector


-- | Calculates the determinant of a 2x2 matrix.
--
--   Pattern matches a 2x2 matrix by its elements and then uses the determinant
--   equation.
det2x2 :: RealFloat a => Matrix a -> a
det2x2 [[a, b], [c, d]] = a * d - b * c
det2x2 _ = error "matrix must be a 2x2 matrix"


-- | Finds the inverse matrix for a 2x2 matrix.
--
--   Pattern matches a 2x2 matrix by its elements, then uses its determinant
--   and adjugate matrix to find the inverse by its mathematical definition.
inv2x2 :: RealFloat a => Matrix a -> Matrix a
inv2x2 matrix@[[a, b], [c, d]]
    | det == 0  = error "matrix is singular"
    | otherwise = mapMatrix (/det) adjugate
    where
        det      = det2x2 matrix
        adjugate = [[d, -b], [-c, a]]
inv2x2 _ = error "matrix must be a 2x2 matrix"


-- | Finds the lower triangular matrix of a Cholesky decomposition for a 2x2
--   matrix.
--
--   Pattern matches a 2x2 matrix by its elements, then uses a decomposed
--   Cholesky–Banachiewicz algorithm for just a 2x2 matrix.
cholesky2x2 :: RealFloat a => Matrix a -> Matrix a
cholesky2x2 [[a, b], [c, d]] = cholesky
    where
        l11      = sqrt a
        l21      = c / l11
        l22      = sqrt (d - l21^^2)
        cholesky = [[l11, 0], [l21, l22]]
cholesky2x2 _ = error "matrix must be a 2x2 matrix"


-- | Performs the outer product of a vector with itself. (x * x^T)
--
--   For each value in the vector, multiply over the entire vector to create a
--   matrix.
outerProd :: RealFloat a => Vector a -> Matrix a
outerProd vector = [map (*v) vector | v <- vector]


-- | Computes the quadratic form of a matrix and vector.  (x^T * A * x)
--
--   Computes (y = A * x) then (x^T * y) by doing matrix multiplication, then
--   the dot product with that result.
quadForm :: RealFloat a => Matrix a -> Vector a -> a
quadForm matrix vector
    | nRows matrix /= nCols matrix
        = error "matrix must be square"
    | nRows matrix /= length vector
        = error "vector must be the same length as rows/columns in matrix"
    | otherwise
        = transVecMatrixVec
    where
        matrixVec         = map sum $ zipRowsWith (*) matrix vector
        transVecMatrixVec = sum $ zipWith (*) vector matrixVec
