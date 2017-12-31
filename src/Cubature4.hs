module Cubature4
  where
import Internal4
import Simplex
import qualified Data.Vector.Unboxed         as UV
import Control.Monad.ST
import Data.Array.Unboxed
import Data.Array.ST as STA

simplicesToArray3 :: Simplices -> ST s (STUArray s (Int,Int,Int) Double) -- U3dArray
simplicesToArray3 simplices = do
  let dim = length (head (head simplices))
      nsimplices = length simplices
      assocList = map (\[i,j,k] -> ((i,j,k), (simplices!!(k-1))!!(j-1)!!(i-1)))
                      (sequence [[1..dim], [1..(dim+1)], [1..nsimplices]])
      arr = array ((1,1,1),(dim,dim+1,nsimplices)) assocList :: UArray (Int,Int,Int) Double
  STA.thaw arr

integrateOnSimplex
    :: (UVectorD -> UVectorD) -- integrand
    -> Simplices              -- domain
    -> Int                    -- number of components
    -> Int                    -- maximum number of evaluations
    -> Double                 -- desired absolute error
    -> Double                 -- desired relative error
    -> Int                    -- integration rule: 1, 2, 3 or 4
    -> ST s ([Double], [Double], Int, Bool) -- integral, error, nevals, code
integrateOnSimplex f s ncomp maxevals absError relError rule = do
  let n = length (head s) - 1
  v <- simplicesToArray3 s
  (vals, errors, nevals, fl) <-
    adsimp n ncomp maxevals f absError relError rule v False
  return (UV.toList vals, UV.toList errors, nevals, fl)

fExample :: UVectorD -> UVectorD
fExample v = let list = UV.toList v in UV.fromList [sum list, sum (map (^2) list)]

fExample' :: UVectorD -> UVectorD
fExample' v = let list = UV.toList v in UV.fromList [sum list, sum list]

example rule = integrateOnSimplex fExample [canonicalSimplex 3] 2 10000 0 1e-5 rule

example' rule = integrateOnSimplex fExample' [canonicalSimplex 3] 2 10000 0 1e-5 rule

fExample2 :: UVectorD -> UVectorD
fExample2 v = UV.singleton $ sqrt((x!!3-x!!2)/(x!!1-x!!0))*exp(-(x!!1-x!!0))
  where y = UV.toList v
        x = map (\i -> sum $ take i y) [1..4]

example2 maxevals rule = integrateOnSimplex fExample2 [canonicalSimplex 4] 1 maxevals 0 1e-5 rule

fExample2' :: UVectorD -> UVectorD
fExample2' x = UV.singleton $ sqrt((x UV.!3 - x UV.!2)/(x UV.!1 - x UV.!0))*exp(-(x UV.!1 - x UV.!0))

example2' maxevals rule = integrateOnSimplex fExample2'
                          [[[0,0,0,0],[1,1,1,1],[0,1,1,1],[0,0,1,1],[0,0,0,1]]]
                          1 maxevals 0 1e-5 rule
