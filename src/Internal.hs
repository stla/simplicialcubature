module Internal
   where
-- import qualified Data.Array                  as A
--import qualified Data.Array.IArray           as IA
import           Data.Array.IArray           (Array)
import           Data.Array.IO               (IOUArray, getBounds, getElems,
                                              newArray, newArray_, newListArray,
                                              readArray, writeArray, mapIndices)
import qualified Data.Array.IO               as IOA
--import qualified Data.Array.MArray           as MA
--import qualified Data.Array.Unboxed          as AU
import           Data.Vector.Unboxed         (Vector, freeze, fromList, toList)
import qualified Data.Vector.Unboxed         as UV
import           Data.Vector.Unboxed.Mutable (IOVector, write, new)
import qualified Data.Vector.Unboxed.Mutable as UMV
import           Math.Combinat.Permutations  (permuteMultiset)
import Control.Monad ((=<<))

type IOMatrix = IOUArray (Int,Int) Double
type IO1dArray = IOUArray Int Double
type Matrix = Array (Int,Int) Double
type IOVectorD = IOVector Double
type IOVectorI = IOVector Int
type UVectorD = Vector Double
type UVectorI = Vector Int

smprms :: Int -> Int -> IO (IOMatrix, IOMatrix, IOVectorI)
smprms n key = do
  let (rls, gms, wts) | key == 1 = (3, 2, 3) :: (Int, Int, Int)
                      | key == 2 = (5, 4, 6) :: (Int, Int, Int)
                      | key == 3 = (7, 7, 11) :: (Int, Int, Int)
                      | key == 4 = (if n == 2 then (7, 11, 20) else (7, 12, 21))
                                    :: (Int, Int, Int)
  w <- newArray ((1,1), (wts, rls)) 0 :: IO IOMatrix
  pts <- UMV.replicate wts 0 :: IO IOVectorI
  g <- newArray ((1,1), (n+1, wts)) 0 :: IO IOMatrix
  let np = n+1
      n2 = np * (n+2)
      n4 = n2 * (n+3) * (n+4)
      n6 = n4 * (n+5) * (n+6)
      n8 = n6 * (n+7) * (n+8)
      ndouble = fromIntegral n :: Double
      r1 = (ndouble + 4 - sqrt(15)) / (ndouble*ndouble + 8*ndouble + 1)
      s1 = 1 - ndouble*r1
      l1 = s1 - r1
  _ <- mapM (\j -> writeArray g (j,1) (1/(ndouble+1))) [1 .. np]
  write pts 0 1
  writeArray g (1,gms+1) s1
  _ <- mapM (\j -> writeArray g (j,gms+1) r1) [2 .. np]
  write pts gms np
  case key < 4 of
    True -> do
      writeArray w (1,rls) 1
      writeArray w (gms+1,rls-1) (1/(ndouble+1))
    False -> return ()
  let iw = if key < 4 then rls-2 else rls
  writeArray g (1,2) (3/(ndouble+3))
  _ <- mapM (\j -> writeArray g (j,2) (1/(ndouble+3))) [2 .. np]
  write pts 1 np
  let n2double = fromIntegral n2 :: Double
  writeArray w (2,iw) ((ndouble+3)^3 / (4*n2double*(ndouble+3)))
  case key > 1 of
    True -> do
      case n == 2 of
        True -> do
          let l2 = 0.62054648267200632589046034361711
              l1 = -sqrt(0.5 - l2*l2)
              r1 = (1-l1)/3
              s1 = 1 - 2*r1
          write pts gms 3
          writeArray w (gms+1,iw-1) (1/6) -- rq: déjà écrit
          let r2 = (1-l2)/3
              s2 = 1 - 2*r2
          writeArray g (1,gms+2) s2
          _ <- mapM (\j -> writeArray g (j,gms+2) r2) [2 .. np]
          write pts (gms+1) 3
          writeArray w (gms+2,iw-1) (1/6)
        False -> do
          let r2 = (ndouble+4+sqrt(15)) / (ndouble*ndouble+8*ndouble+1)
              s2 = 1 - ndouble*r2
              l2 = s2 - r2
          writeArray g (1,gms+2) s2
          _ <- mapM (\j -> writeArray g (j,gms+2) r2) [2 .. np]
          write pts (gms+1) np
          writeArray w (gms+2,iw-1) ((2/(ndouble+3)-l1)/(n2double*(l2-l1)*l2*l2))
          writeArray w (gms+1,iw-1) ((2/(ndouble+3)-l2)/(n2double*(l1-l2)*l1*l1))
      writeArray g (1,3) (5/(ndouble+5))
      _ <- mapM (\j -> writeArray g (j,3) (1/(ndouble+5))) [2 .. np]
      write pts 2 np
      writeArray g (1,4) (3/(ndouble+5))
      writeArray g (2,4) (3/(ndouble+5))
      _ <- mapM (\j -> writeArray g (j,4) (1/(ndouble+5))) [3 .. np]
      write pts 3 (div (n*np) 2)
      let n4double = fromIntegral n4 :: Double
      writeArray w (2,iw-2) (-(ndouble+3)^5 / (16*n4double))
      writeArray w (3,iw-2) ((ndouble+5)^5 / (16*n4double*(ndouble+5)))
      writeArray w (4,iw-2) ((ndouble+5)^5 / (16*n4double*(ndouble+5)))
    False -> return ()
  case key > 2 of
    True -> do
      let u1 = (ndouble+7+2*sqrt(15)) / (ndouble^2+14*ndouble-11)
          v1 = u1 --- XXX
          d1 = v1 - u1
      writeArray g (1,gms+3) v1
  return (g, w, pts)

matprod :: IOMatrix -> [Double] -> IO UVectorD
matprod mat x = do
  (_, (m,n)) <- getBounds mat
  out <- UMV.new m :: IO IOVectorD
  let step i | i == m+1 = freeze out
             | otherwise = do
              coef <- innerstep i 1 0
              write out (i-1) coef
              step (i+1)
      innerstep :: Int -> Int -> Double -> IO Double
      innerstep i j s | j == n+1 = return s
                      | otherwise = do
                       mat_ij <- readArray mat (i,j)
                       innerstep i (j+1) (s + mat_ij * (x !! (j-1)))
  step 1

smpsms :: IOMatrix -> Int -> (UVectorD -> UVectorD) -> IO1dArray
       -> Double -> IO IO1dArray
smpsms vertex nf f g scalar = do
  gAsList <- getElems g
  f_gPermutations <- mapM ((fmap toList) . (fmap f) . (matprod vertex))
                          (permuteMultiset gAsList)
  newListArray (1,nf) (map ((*) scalar) (foldl1 (zipWith (+)) f_gPermutations))

extractColumn :: IOMatrix -> Int -> IO IO1dArray
extractColumn m j = do
  (_, (nrow,_)) <- getBounds m
  mapIndices (1,nrow) (\i -> (i,j)) m

extractRow :: IOMatrix -> Int -> IO IO1dArray
extractRow m i = do
  (_, (_,ncol)) <- getBounds m
  mapIndices (1,ncol) (\j -> (i,j)) m

outerProduct :: IO1dArray -> IO1dArray -> IO IOMatrix
outerProduct x1 x2 = do
  (_, n1) <- getBounds x1
  (_, n2) <- getBounds x2
  out <- newArray_ ((1,1),(n1,n2)) :: IO IOMatrix
  let step :: Int -> Int -> IO IOMatrix
      step i j | i == n1 && j == n2+1 = return out
               | j == n2+1 = step (i+1) 1
               | otherwise = do
                  x1_i <- readArray x1 i
                  x2_j <- readArray x2 j
                  writeArray out (i,j) (x1_i * x2_j)
                  step i (j+1)
  step 1 1

testouterProduct :: IO Matrix
testouterProduct = do
  x1 <- newListArray (1,2) [1, 2] :: IO IO1dArray
  x2 <- newListArray (1,3) [1, 2, 3] :: IO IO1dArray
  (>>=) (outerProduct x1 x2) IOA.freeze

sumMatrices :: [IOMatrix] -> IO IOMatrix
sumMatrices matrices = do
  (_, (n1,n2)) <- getBounds (head matrices)
  out <- newArray_ ((1,1),(n1,n2)) :: IO IOMatrix
  let step :: Int -> Int -> IO IOMatrix
      step i j | i == n1 && j == n2+1 = return out
               | j == n2+1 = step (i+1) 1
               | otherwise = do
                 coefs <- mapM (\m -> readArray m (i,j)) matrices
                 writeArray out (i,j) (sum coefs)
                 step i (j+1)
  step 1 1

smprul :: IOMatrix -> Int -> (UVectorD -> UVectorD) -> Double -> IOMatrix
       -> IOMatrix -> UVectorI -> IO (IO1dArray, IO1dArray)
smprul vrts nf f vol g w pts = do
  let rtmn = 0.1
      small = 1e-12
      errcof = 8
  (_, (_,rls)) <- getBounds w
  let ptsPositive = toList $ UV.map (+1) $ UV.findIndices (> 0) pts
  toSum <- mapM (\k -> do
                         g_colk <- extractColumn g k
                         w_rowk <- extractRow w k
                         sms <- smpsms vrts nf f g_colk vol
                         outerProduct sms w_rowk)
                ptsPositive
  rule <- sumMatrices toSum -- quid si la liste est vide ?
  basval <- extractColumn rule 1
  rgnerr <- newArray (1,nf) 0 :: IO IO1dArray
  let step :: Int -> IO ()
      step i | i == nf+1 = return ()
             | otherwise = do
                basval_i <- readArray basval i
                let nmbs = abs basval_i
                (rt, nmcp) <- innerstep rls rtmn nmbs
                case rt < 1 && rls > 3 of
                  True  -> writeArray rgnerr i (rt*nmcp)
                  False -> return ()
                rgnerr_i <- readArray rgnerr i
                writeArray rgnerr i (max (errcof*rgnerr_i) (small*nmbs))
                step (i+1)
              where
               innerstep :: Int -> Double -> Double -> IO (Double, Double)
               innerstep k x y | k == 1 = return (x, y)
                               | otherwise = do
                                rule_ik <- readArray rule (i,k)
                                rule_ikm1 <- readArray rule (i,k-1)
                                let nmrl = max (abs rule_ik) (abs rule_ikm1)
                                rgnerr_i <- readArray rgnerr i
                                writeArray rgnerr i (max nmrl rgnerr_i)
                                case nmrl > small*y && k < rls of
                                  True -> innerstep (k-2) (max (nmrl/y) x) nmrl
                                  False -> innerstep (k-2) x nmrl
  step 1
  return (basval, rgnerr)

testsmprul :: IO ([Double], [Double])
testsmprul = do
  let nd = 3
  vrts <- newListArray ((1,1),(nd,nd+1)) [1,2,3,4,5,6,7,8,9,10,11,12] :: IO IOMatrix
  let nf = 2
  let f :: UVectorD -> UVectorD
      f v = let list = toList v in fromList [sum list, product list]
  let vol = 2.0
  g <- newListArray ((1,1),(nd+1,3)) [1,2,3,4,5,6,7,8,9,10,11,12] :: IO IOMatrix
  w <- newListArray ((1,1),(3,3)) [1,2,3,4,5,6,7,8,9] :: IO IOMatrix
  let pts = fromList [1,2,3] :: UVectorI
  -- gcol <- extractColumn g 1
  -- sms <- smpsms vrts nf f gcol 1
  -- smsAsList <- getElems sms
  -- return (smsAsList, smsAsList)
  (basval, rgnerr) <- smprul vrts nf f vol g w pts
  basvalList <- getElems basval
  rgnerrList <- getElems rgnerr
  return (basvalList, rgnerrList)

rowMeans :: IOMatrix -> IO UVectorD
rowMeans m = do
  (_, (nrow,ncol)) <- getBounds m
  rowSumsIO <- new nrow :: IO IOVectorD
  let step :: Int -> IO ()
      step i | i == nrow+1 = return ()
             | otherwise = do
                sum_i <- inner 1 0
                write rowSumsIO (i-1) sum_i
                step (i+1)
              where
                inner :: Int -> Double -> IO Double
                inner j x | j == ncol+1 = return x
                          | otherwise = do
                             coef <- readArray m (i,j)
                             inner (j+1) (x + coef)
  step 1
  rowSums <- freeze rowSumsIO
  return $ UV.map (/(fromIntegral ncol)) rowSums

testrowMeans :: IO UVectorD
testrowMeans = do
  m <- newListArray ((1,1),(2,3)) [1,2,3,4,5,6] :: IO IOMatrix
  rowMeans m

array1dToUVectorD :: IO1dArray -> IO UVectorD
array1dToUVectorD array = do
  (=<<) (return . fromList) (getElems array)

smpdfs :: Int -> Int -> (UVectorD -> UVectorD) -> Int -> Int
       -> IOUArray (Int,Int,Int) Double -> IO (Int, IOUArray (Int,Int,Int) Double)
smpdfs nd nf f top sbs vrts = do
  let cuttf = 2.0
      cuttb = 8.0
      --is = 1
      --js = 2
      --dfmx = 0.0
--      emx = 0.0
  v <- mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,top)) vrts
  cn <- rowMeans v
  let fc = f cn
      dfmd = sum $ toList $ UV.map abs fc
  frthdf <- newArray ((1,1),(nd,nd+1)) 0 :: IO IOMatrix
  iejeitjtisjs <- new 6 :: IO IOVectorI
  write iejeitjtisjs 4 1
  write iejeitjtisjs 5 2
  dfmxdfnx <- UMV.replicate 2 0 :: IO IOVectorD
  let step :: Int -> Double -> IO ()
      step i x | i == nd+1 = return ()
               | otherwise = do
                  emx <- inner (i+1) x
                  step (i+1) emx
                where
                  inner :: Int -> Double -> IO Double
                  inner j y | j == nd+2 = return y
                            | otherwise = do
                              viIO <- extractColumn v i
                              vjIO <- extractColumn v j
                              vi <- array1dToUVectorD viIO
                              vj <- array1dToUVectorD vjIO
                              let h = UV.map (*(2/(5*(fromIntegral nd +1))))
                                             (UV.zipWith (-) vi vj)
                                  ewd = sum $ toList $ UV.map abs h
                                  twoh = UV.map (*2) h
                                  t1 = f (UV.zipWith (-) cn twoh)
                                  t2 = f (UV.zipWith (+) cn twoh)
                                  t3 = UV.map (*6) fc
                                  t4 = f (UV.zipWith (-) cn h)
                                  t5 = f (UV.zipWith (+) cn h)
                                  t6 = UV.map (*(-4)) (UV.zipWith (+) t4 t5)
                                  tsum = (foldl1 (UV.zipWith (+))) [t1,t2,t3,t6]
                                  dfr1 = sum $ toList $ UV.map abs tsum
                                  dfr2 = if dfmd+dfr1/8 == dfmd then 0 else dfr1
                                  dfr3 = dfr2*ewd
                              dfmx <- UMV.read dfmxdfnx 0
                              case dfr3 >= dfmx of
                                True -> do
                                  is <- UMV.read iejeitjtisjs 4
                                  js <- UMV.read iejeitjtisjs 5
                                  write iejeitjtisjs 2 is
                                  write iejeitjtisjs 3 js
                                  write iejeitjtisjs 4 i
                                  write iejeitjtisjs 5 j
                                  write dfmxdfnx 1 dfmx
                                  write dfmxdfnx 0 dfr3
                                False -> do
                                  dfnx <- UMV.read dfmxdfnx 1
                                  case dfr3 >= dfnx of
                                    True -> do
                                      write iejeitjtisjs 2 i
                                      write iejeitjtisjs 3 j
                                      write dfmxdfnx 1 dfr3
                                    False -> return ()
                              writeArray frthdf (i,j) dfr3
                              case ewd >= y of
                                True -> do
                                  write iejeitjtisjs 0 i
                                  write iejeitjtisjs 1 j
                                  inner (j+1) ewd
                                False -> inner (j+1) y
  step 1 0
  dfmx <- UMV.read dfmxdfnx 0
  dfnx <- UMV.read dfmxdfnx 1
  let nregions = if dfnx > dfmx/cuttf then 4 else 3
  case dfnx > dfmx/cuttf of
    True -> return ()
    False -> do
      case dfmx == 0 of
        True -> do
          ie <- UMV.read iejeitjtisjs 0
          je <- UMV.read iejeitjtisjs 1
          write iejeitjtisjs 4 ie
          write iejeitjtisjs 5 je
        False -> do
          let loop :: Int -> Double -> Int -> IO Int
              loop l x ls | l == nd+2 = return ls
                          | otherwise = do
                            is <- UMV.read iejeitjtisjs 4
                            js <- UMV.read iejeitjtisjs 5
                            case (l /= is) && (l /= js) of
                              True -> do
                                let it = minimum [l,is,js]
                                    jt = maximum [l,is,js]
                                write iejeitjtisjs 2 it
                                write iejeitjtisjs 3 jt
                                let lt = is+js+l-it-jt
                                dfr1 <- readArray frthdf (it,lt)
                                dfr2 <- readArray frthdf (lt,jt)
                                let dfr = dfr1 + dfr2
                                case dfr >= x of
                                  True -> loop (l+1) dfr l
                                  False -> loop (l+1) x ls
                              False -> loop (l+1) x ls
          ls <- loop 1 0 0
          is <- UMV.read iejeitjtisjs 4
          js <- UMV.read iejeitjtisjs 5
          difil <- readArray frthdf (min is ls, max is ls)
          diflj <- readArray frthdf (min js ls, max js ls)
          let dfnx = max difil diflj
          write dfmxdfnx 1 dfnx
          case dfmx/cuttb < dfnx && difil > diflj of
            True -> do
              is <- UMV.read iejeitjtisjs 4
              js <- UMV.read iejeitjtisjs 5
              it <- UMV.read iejeitjtisjs 2
              write iejeitjtisjs 2 is
              write iejeitjtisjs 4 js
              write iejeitjtisjs 5 it
            False -> return ()
  vrts2 <- newArray_ ((1,1,1),(nd,nd+1,sbs+nregions-1)) :: IO (IOUArray (Int,Int,Int) Double)
  let go1 :: Int -> IO ()
      go1 i | i == nd+1 = return ()
            | otherwise = do
              let go2 j | j == nd+2 = go1 (i+1)
                        | otherwise = do
                          let go3 k | k == sbs+nregions = go2 (j+1)
                                    | otherwise = do
                                      case k <= sbs of
                                        True -> do
                                          coef <- readArray vrts (i,j,k)
                                          writeArray vrts2 (i,j,k) coef
                                        False -> do
                                          coef <- readArray v (i,j)
                                          writeArray vrts2 (i,j,k) coef
                                      go3 (k+1)
                          go3 1
              go2 1
  go1 1
  is <- UMV.read iejeitjtisjs 4
  js <- UMV.read iejeitjtisjs 5
  vtiIO <- extractColumn v is
  vtjIO <- extractColumn v js
  vti <- array1dToUVectorD vtiIO
  vtj <- array1dToUVectorD vtjIO
  case nregions == 4 of
    True -> do
      let vt = UV.map (/2) (UV.zipWith (+) vti vtj)
      replaceDimension vrts2 (js,top) vt
      replaceDimension vrts2 (is,sbs+1) vti
    False -> do
      let vt = UV.map (/3) (UV.zipWith (+) (UV.map (*2) vti) vtj)
      return ()
  return (nregions, vrts2)

replaceDimension :: IOUArray (Int,Int,Int) Double -> (Int,Int) -> UVectorD -> IO ()
replaceDimension m (j,k) v = do
  (_, (n,_,_)) <- getBounds m
  let loop :: Int -> IO ()
      loop i | i == n+1 = return ()
             | otherwise = do
               writeArray m (i,j,k) ((UV.!) v (i-1))
               loop (i+1)
  loop 1
