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
-- ok si je réduis maxEvals :
-- *Cubature> integrateOnSimplex fExample [canonicalSimplex 3] 2 10000 0 1e-5 1
-- ([1.485486158822895,0.48644692253316674],[1.485486158822925e-12,0.2675948891062603],9954,True)
-- *Cubature> integrateOnSimplex fExample [canonicalSimplex 3] 2 30 0 1e-5 1
-- ([0.125,4.9999999999999996e-2],[1.25e-13,0.4048388923847347],9,True)
-- => instabilité numérique ?
fExample2 :: UVectorD -> UVectorD
fExample2 v = UV.singleton $ sqrt((x!!3-x!!2)/(x!!1-x!!0))*exp(-(x!!1-x!!0))
  where y = UV.toList v
        x = map (\i -> sum $ take i y) [1..4]

example2 maxevals rule = integrateOnSimplex fExample2 [canonicalSimplex 4] 1 maxevals 0 1e-5 rule
