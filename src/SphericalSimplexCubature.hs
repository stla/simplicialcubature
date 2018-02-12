module SphericalSimplexCubature
  where
import           Data.Vector.Unboxed               (Vector, toList)
import           Simplex
import           SimplexCubature
import           SphericalSimplexCubature.Internal

integrand :: [Double] -> Double
integrand x = x!!0 * x!!0 * x!!3 + x!!1 * x!!1 * x!!3 + x!!2 * x!!2 * x!!3 + x!!3 * x!!3 * x!!3

ssimplex :: SphericalSimplex
ssimplex = orthants 4 !! 15

integrand' :: [Double] -> Double
integrand' = transformedIntegrand ssimplex integrand

integrand'' :: Vector Double -> Double
integrand'' = integrand' . toList

test :: IO Result
test = integrateOnSimplex' integrand'' [canonicalSimplex 3] 100000 0 1e-5 3

-- sphere surface
integrand2 :: [Double] -> Double
integrand2 _ = 1

ssimplices :: [SphericalSimplex]
ssimplices = orthants 3

integrands :: [[Double] -> Double]
integrands = map (`transformedIntegrand` integrand2) ssimplices

integrands' :: [Vector Double -> Double]
integrands' = map (. toList) integrands

integrals :: IO [Result]
integrals = mapM (\f -> integrateOnSimplex' f [canonicalSimplex 2] 10000 0 1e-5 3) integrands'

total :: IO Double
total = do
  results <- integrals
  return $ sum (map value results)
