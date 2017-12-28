module Simplex
 where
import           Common
import           Data.Array.IO      (thaw)
import           Data.Array.Unboxed (UArray, array)
import           Data.Matrix        (detLU, elementwiseUnsafe, fromLists,
                                     multStd)

type Simplex = [[Double]]
type Simplices = [Simplex]

isValidSimplex :: Simplex -> Bool
isValidSimplex simplex =
  (length simplex == dim+1) &&
    (all (== dim) (map length (tail simplex)))
  where dim = length (head simplex)

isValidSimplices :: Simplices -> Bool
isValidSimplices simplices =
  (foldr (&&) True (map isValidSimplex simplices)) &&
    (all (== spaceDim (head simplices))) (map spaceDim (tail simplices))
  where spaceDim simplex = length (head simplex)

simplicesToArray :: Simplices -> IO IO3dArray
simplicesToArray simplices = do
  let dim = length (head (head simplices))
      nsimplices = length simplices
      assocList = map (\[i,j,k] -> ((i,j,k), (simplices !! k) !! j !! i))
                      (sequence [[1..dim], [1..(dim+1)], [1..nsimplices]])
      arr = array ((1,1,1),(dim,dim+1,nsimplices)) assocList
            :: UArray (Int,Int,Int) Double
  thaw arr

canonicalSimplex :: Int -> Simplex
canonicalSimplex dim =
  (replicate dim 0) :
    (map (\v -> map (fromIntegral.fromEnum.(== v)) [1..dim]) [1..dim])

simplexVolume :: Simplex -> Double
simplexVolume s = (abs (detLU v)) / (fromIntegral (product [1..n]))
  where n = length s - 1
        m1 = fromLists (tail s)
        m2 = fromLists $ replicate n (head s)
        v = elementwiseUnsafe (-) m1 m2

jacobian :: Simplex -> Double
jacobian s = abs (detLU (elementwiseUnsafe (-) m1 m2))
  where m1 = fromLists (tail s)
        m2 = fromLists $ replicate (length s - 1) (head s)
