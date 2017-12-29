module Cubature
  where
import Internal
import Simplex
import qualified Data.Vector.Unboxed         as UV

integrateOnSimplex
    :: (UVectorD -> UVectorD) -- integrand
    -> Simplices              -- domain
    -> Int                    -- number of components
    -> Int                    -- maximum number of evaluations
    -> Double                 -- desired absolute error
    -> Double                 -- desired relative error
    -> Int                    -- integration rule: 1, 2, 3 or 4
    -> IO ([Double], [Double], Int, Bool) -- integral, error, nevals, code
integrateOnSimplex f s ncomp maxevals absError relError rule = do
  let n = length (head s) - 1
  v <- simplicesToArray s
  (vals, errors, nevals, fl) <-
    adsimp n ncomp maxevals f absError relError rule v False
  return (UV.toList vals, UV.toList errors, nevals, fl)

fExample :: UVectorD -> UVectorD
fExample v = let list = UV.toList v in UV.fromList [sum list, sum (map (^2) list)]

fExample' :: UVectorD -> UVectorD
fExample' v = let list = UV.toList v in UV.fromList [sum list, sum list]

example rule = integrateOnSimplex fExample [canonicalSimplex 3] 2 10000 0 1e-5 rule

example' rule = integrateOnSimplex fExample' [canonicalSimplex 3] 2 10000 0 1e-5 rule

-- rule 2 ok
-- rule 1 bizarre si je change une composante ça change le résultat de l'autre
