module Internal
  where
import           Common
import           Data.Array.IArray           (Array)
import           Data.Array.IO               (IOUArray, getBounds, getElems,
                                              mapIndices, newArray, newArray_,
                                              newListArray, readArray,
                                              writeArray)
import qualified Data.Array.IO               as IOA
--import qualified Data.Array.MArray           as MA
-- import qualified Data.Array                  as A
--import qualified Data.Array.IArray           as IA
import qualified Data.Array.Unboxed          as UA
import           Data.Array.Unboxed          (UArray)
import           Control.Monad               ((=<<))
import           Data.Vector.Unboxed         (Vector, freeze, fromList, toList)
import qualified Data.Vector.Unboxed         as UV
import qualified Data.Vector as V
import           Data.Vector.Unboxed.Mutable (IOVector, new, write)
import qualified Data.Vector.Unboxed.Mutable as UMV
import           Math.Combinat.Permutations  (permuteMultiset)
import           Simplex                     (Simplex, simplexVolume)
import           Data.Maybe                  (fromJust)
import           Data.Sequence               (Seq, update, (><), index)
import qualified Data.Sequence               as S
import qualified Data.Foldable as DF

type IOMatrix = IOUArray (Int,Int) Double
type IO1dArray = IOUArray Int Double
type IMatrix = Array (Int,Int) Double
type UMatrix = UArray (Int,Int) Double
type IOVectorD = IOVector Double
type IOVectorI = IOVector Int
type UVectorD = Vector Double
type UVectorI = Vector Int

fromInt :: Int -> Double
fromInt = fromIntegral

smprms :: Int -> Int -> IO (IOMatrix, [UVectorD], IOVectorI)
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
      ndbl = fromInt n
      r1 = (ndbl + 4 - sqrt(15)) / (ndbl*ndbl + 8*ndbl + 1)
      s1 = 1 - ndbl*r1
      l1 = s1 - r1
  _ <- mapM (\j -> writeArray g (j,1) (1/(ndbl+1))) [1 .. np]
  write pts 0 1
  writeArray g (1,gms+1) s1
  _ <- mapM (\j -> writeArray g (j,gms+1) r1) [2 .. np]
  write pts gms np
  case key < 4 of
    True -> do
      writeArray w (1,rls) 1
      writeArray w (gms+1,rls-1) (1/(ndbl+1))
    False -> return ()
  let iw = if key < 4 then rls-2 else rls
  writeArray g (1,2) (3/(ndbl+3))
  _ <- mapM (\j -> writeArray g (j,2) (1/(ndbl+3))) [2 .. np]
  write pts 1 np
  let n2double = fromInt n2 :: Double
  writeArray w (2,iw) ((ndbl+3)^3 / (4*n2double*(ndbl+3)))
  case key > 1 of
    True -> do
      case n == 2 of
        True -> do
          let l2 = 0.62054648267200632589046034361711
              l1 = -sqrt(0.5 - l2*l2)
              r1 = (1-l1)/3
              s1 = 1 - 2*r1
          writeArray g (1,gms+1) s1
          _ <- mapM (\j -> writeArray g (j,gms+1) r1) [2 .. np]
          write pts gms 3
          writeArray w (gms+1,iw-1) (1/6) -- rq: déjà écrit
          let r2 = (1-l2)/3
              s2 = 1 - 2*r2
          writeArray g (1,gms+2) s2
          _ <- mapM (\j -> writeArray g (j,gms+2) r2) [2 .. np]
          write pts (gms+1) 3
          writeArray w (gms+2,iw-1) (1/6)
        False -> do
          let r2 = (ndbl+4+sqrt(15)) / (ndbl*ndbl+8*ndbl+1)
              s2 = 1 - ndbl*r2
              l2 = s2 - r2
          writeArray g (1,gms+2) s2
          _ <- mapM (\j -> writeArray g (j,gms+2) r2) [2 .. np]
          write pts (gms+1) np
          writeArray w (gms+2,iw-1) ((2/(ndbl+3)-l1)/(n2double*(l2-l1)*l2*l2))
          writeArray w (gms+1,iw-1) ((2/(ndbl+3)-l2)/(n2double*(l1-l2)*l1*l1))
      writeArray g (1,3) (5/(ndbl+5))
      _ <- mapM (\j -> writeArray g (j,3) (1/(ndbl+5))) [2 .. np]
      write pts 2 np
      writeArray g (1,4) (3/(ndbl+5))
      writeArray g (2,4) (3/(ndbl+5))
      _ <- mapM (\j -> writeArray g (j,4) (1/(ndbl+5))) [3 .. np]
      write pts 3 (div (n*np) 2)
      let n4double = fromInt n4 :: Double
      writeArray w (2,iw-2) (-(ndbl+3)^5 / (16*n4double))
      writeArray w (3,iw-2) ((ndbl+5)^5 / (16*n4double*(ndbl+5)))
      writeArray w (4,iw-2) ((ndbl+5)^5 / (16*n4double*(ndbl+5)))
    False -> return ()
  case key > 2 of
    True -> do
      let u1 = (ndbl+7+2*sqrt(15)) / (ndbl^2+14*ndbl-11)
          v1 = (1-(ndbl-1)*u1)/2
          d1 = v1 - u1
      writeArray g (1,gms+3) v1
      writeArray g (2,gms+3) v1
      _ <- mapM (\j -> writeArray g (j,gms+3) u1) [3 .. np]
      write pts (gms+2) (div (n*np) 2)
      let u2 = (ndbl+7-2*sqrt(15)) / (ndbl^2+14*ndbl-11)
          v2 = (1-(ndbl-1)*u2)/2
          d2 = v2 - u2
      writeArray g (1,gms+4) v2
      writeArray g (2,gms+4) v2
      _ <- mapM (\j -> writeArray g (j,gms+4) u2) [3 .. np]
      write pts (gms+3) (div (n*np) 2)
      case n == 2 of
        True -> do
          writeArray w (gms+3,iw-3) ((155-sqrt(15))/1200)
          writeArray w (gms+4,iw-3) ((155+sqrt(15))/1200)
          writeArray w (1,iw-3) (1 - 31/40)
        False -> do
          case n == 3 of
            True -> do
              writeArray w (gms+1,iw-3) ((2665+14*sqrt(15))/37800)
              writeArray w (gms+2,iw-3) ((2665-14*sqrt(15))/37800)
              writeArray w (gms+3,iw-3) (2*15/567)
              write pts (gms+3) 0
            False -> do
              let r2 = (ndbl+4+sqrt(15)) / (ndbl*ndbl+8*ndbl+1)
                  l2 = 1 - (ndbl+1)*r2
                  den = l1^4*(l1-l2)*(fromInt n4)
              writeArray w (gms+1,iw-3)
                           ((2*(27-ndbl)/(ndbl+5)-l2*(13-ndbl))/den)
              writeArray w (gms+2,iw-3)
                           ((2*(27-ndbl)/(ndbl+5)-l1*(13-ndbl))/den)
              writeArray w (gms+3,iw-3)
                           ((2/(ndbl+5)-d2)/((fromInt n4)*(d1-d2)*d1^4))
              writeArray w (gms+4,iw-3)
                           ((2/(ndbl+5)-d1)/((fromInt n4)*(d2-d1)*d2^4))
      writeArray g (1,5) (7/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,5) (1/(ndbl+7))) [2 .. np]
      write pts 4 np
      writeArray g (1,6) (5/(ndbl+7))
      writeArray g (2,6) (3/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,6) (1/(ndbl+7))) [3 .. np]
      write pts 5 (np*n)
      _ <- mapM (\i -> writeArray g (i,7) (3/(ndbl+7))) [1,2,3]
      write pts 6 (div ((n-1)*n*np) 6)
      writeArray w (2,iw-4) ((ndbl+3)^7/(fromInt $ 128*n4*(n+5)))
      _ <- mapM
           (\i -> writeArray w (i,iw-4) ((-(ndbl+5)^7)/(fromInt $ 64*n6))) [3,4]
      _ <- mapM
           (\i -> writeArray w (i,iw-4) (((ndbl+7)^7)/(fromInt $ 64*n6*(n+7))))
           [5,6,7]
      return ()
    False -> return ()
  case key == 4 of
    True -> do
      let sg = 1/(fromInt $ 23328*n6)
          u5 = -6^3 * sg * (fromInt $ 52212 - n*(6353 + n*(1934-n*27)))
          u6 =  6^4 * sg * (fromInt $ 7884 - n*(1541 - n*9))
          u7 = -6^5 * sg * (fromInt $ 8292 - n*(1139 - n*3))/(ndbl + 7)
          p0 = -144 * (142528 + n*(23073 - n*115))
          p1 = -12 * (6690556 + n*(2641189 + n*(245378 - n*1495)))
          p2 = -16 * (6503401 + n*(4020794+n*(787281+n*(47323-n*385))))
          p3 = -(6386660 + n*(4411997+n*(951821+n*(61659-n*665))))*(n + 7)
          a = (fromInt p2)/(fromInt $ 3*p3)
          p = a*((fromInt p1)/(fromInt p2) - a)
          q = a*(2*a*a - (fromInt p1)/(fromInt p3)) + (fromInt p0)/(fromInt p3)
          th = acos(-q/(2*(sqrt(-p*p*p))))/3
          r = 2*sqrt(-p*p*p)**(1/3)
          tp = 2*pi/3
          a1 = -a + r*(cos(th))
          a2 = -a + r*(cos(th+2*tp))
          a3 = -a + r*(cos(th+tp))
          npdbl = fromInt np
      writeArray g (1,gms+5) ((1-ndbl*a1)/npdbl)
      _ <- mapM (\i -> writeArray g (i,gms+5) ((1+a1)/npdbl)) [2..np]
      write pts (gms+4) np
      writeArray g (1,gms+6) ((1-ndbl*a2)/npdbl)
      _ <- mapM (\i -> writeArray g (i,gms+6) ((1+a2)/npdbl)) [2..np]
      write pts (gms+5) np
      writeArray g (1,gms+7) ((1-ndbl*a3)/npdbl)
      _ <- mapM (\i -> writeArray g (i,gms+7) ((1+a3)/npdbl)) [2..np]
      write pts (gms+6) np
      writeArray w (gms+5,iw-5)
                 ((u7-(a2+a3)*u6+a2*a3*u5)/(a1^2-(a2+a3)*a1+a2*a3)/a1^5)
      writeArray w (gms+6,iw-5)
                 ((u7-(a1+a3)*u6+a1*a3*u5)/(a2^2-(a1+a3)*a2+a1*a3)/a2^5)
      writeArray w (gms+7,iw-5)
                 ((u7-(a2+a1)*u6+a2*a1*u5)/(a3^2-(a2+a1)*a3+a2*a1)/a3^5)
      writeArray g (1,gms+8) (4/(ndbl+7))
      writeArray g (2,gms+8) (4/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,gms+8) (1/(ndbl+7))) [3..np]
      write pts (gms+7) (div (np*n) 2)
      writeArray w (gms+8,iw-5) (10*(ndbl+7)^6/(fromInt $ 729*n6))
      writeArray g (1,gms+9) (11/(npdbl+7)/2)
      writeArray g (2,gms+9) (5/(npdbl+7)/2)
      _ <- mapM (\i -> writeArray g (i,gms+9) (1/(ndbl+7))) [3..np]
      write pts (gms+7) (np*n)
      writeArray w (gms+9,iw-5) (64*(ndbl+7)^6/(fromInt $ 6561*n6))
      writeArray w (4,iw-5) ((-(ndbl+5)^7)/(fromInt $ 64*n6))
      writeArray w (7,iw-5) (((ndbl+7)^7)/(fromInt $ 64*n6*(n+7)))
      -- XXX
    False -> return ()
  rowsIO <- mapM (extractRow w) [2..wts]
  rows <- mapM array1dToUVectorD rowsIO
  --let cols = map (\j -> fromList $ map (\i -> rows!!(i-2) UV.!(j-1)) [2..wts]) [1..rls]
  let cols = transpose rows
  ptsU <- (=<<) (return.(UV.map fromInt)) (UV.freeze pts)
  let row1 = map (\j -> 1 - UV.foldr (+) 0 (UV.zipWith (*) (UV.tail ptsU) (cols!!(j-1)))) [1..rls]
--  let wmat = fromLists (row1 : (map toList rows))
      wcols = map (\j -> UV.cons (row1!!j) (cols!!j)) [0..(rls-1)]
      wcols2 = head wcols : (map (\col -> UV.zipWith (-) col (head wcols)) (tail wcols))
      nb = UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map (^2) (head wcols)))
      ratio = nb / (UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map (^2) (wcols2!!1))))
      wcol2 = UV.map (*(sqrt ratio)) (wcols2!!1)
      wcols3 = (head wcols) : (wcol2 : ((tail.tail) wcols2))
      -- goal : foldr updateW wcols3 [3..rls]
  let updateW :: [UVectorD] -> Int -> [UVectorD]
      updateW cols k = [cols!!j | j <- [0..(k-2)]] ++ [wknew] ++ [cols!!j | j <- [k..(rls-1)]]
                       where
                        ptsW = UV.map (/nb) (UV.zipWith (*) ptsU (cols!!(k-1)))
                        slice = [cols!!j | j <- [1..(k-2)]]
                        prod1 = fromList $ map ((UV.foldr (+) 0).(UV.zipWith (*) ptsW)) slice
                        --rows = map (\i -> fromList $ map (\j -> slice!!(j-1) UV.!(i-1)) [1..(k-2)]) [1..wts]
                        rows = transpose slice
                        prod2 = fromList $ map ((UV.foldr (+) 0).(UV.zipWith (*) prod1)) rows
                        wk = UV.zipWith (-) (cols!!(k-1)) prod2
                        ratio = nb / (UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map (^2) wk)))
                        wknew = UV.map (*(sqrt ratio)) wk
      wcolsnew = foldl updateW wcols3 [3..rls]
  return (g, transpose wcolsnew, pts)

transpose :: [UVectorD] -> [UVectorD] -- à utiliser 2 fois avant
transpose cols =
  map (\i -> fromList $ map (\j -> cols!!j UV.!i)
      [0..(length cols - 1)]) [0..(UV.length (head cols) - 1)]


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
  f_gPermuts <- mapM ((fmap f).(matprod vertex)) (permuteMultiset gAsList)
  newListArray (1,nf)
               (toList (UV.map (*scalar) (foldl1 (UV.zipWith (+)) f_gPermuts)))

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

testouterProduct :: IO IMatrix
testouterProduct = do
  x1 <- newListArray (1,2) [1, 2] :: IO IO1dArray
  x2 <- newListArray (1,3) [1, 2, 3] :: IO IO1dArray
  (>>=) (outerProduct x1 x2) IOA.freeze

outerProduct2 :: IO1dArray -> UVectorD -> IO IOMatrix
outerProduct2 x1 x2 = do
  (_, n1) <- getBounds x1
  let n2 = UV.length x2
  out <- newArray_ ((1,1),(n1,n2)) :: IO IOMatrix
  let step :: Int -> Int -> IO IOMatrix
      step i j | i == n1 && j == n2+1 = return out
               | j == n2+1 = step (i+1) 1
               | otherwise = do
                  x1_i <- readArray x1 i
                  writeArray out (i,j) (x1_i * (x2 UV.! (j-1)))
                  step i (j+1)
  step 1 1

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
       -> [UVectorD] -> UVectorI -> IO (IO1dArray, IO1dArray)
smprul vrts nf f vol g w pts = do
  let rtmn = 0.1
      small = 1e-12
      errcof = 8
  --(_, (_,rls)) <- getBounds w
      rls = UV.length (head w)
      ptsPositive = toList $ UV.map (+1) $ UV.findIndices (> 0) pts
  toSum <- mapM (\k -> do
                         g_colk <- extractColumn g k
--                         w_rowk <- extractRow w k
                         sms <- smpsms vrts nf f g_colk vol
                         outerProduct2 sms (w!!(k-1)))
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
  let w = [fromList [1,2,3], fromList [4,5,6], fromList [7,8,9]] :: [UVectorD]
  let pts = fromList [1,2,3] :: UVectorI
  -- gcol <- extractColumn g 1
  -- sms <- smpsms vrts nf f gcol 1
  -- smsAsList <- getElems sms
  -- return (smsAsList, smsAsList)
  (basval, rgnerr) <- smprul vrts nf f vol g w pts
  basvalList <- getElems basval
  rgnerrList <- getElems rgnerr
  return (basvalList, rgnerrList)

rowSums :: IOMatrix -> IO UVectorD
rowSums m = do
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
  freeze rowSumsIO

rowMeans :: IOMatrix -> IO UVectorD
rowMeans m = do
  (_, (_,ncol)) <- getBounds m
  rowsums <- rowSums m
  return $ UV.map (/(fromInt ncol)) rowsums

testrowMeans :: IO UVectorD
testrowMeans = do
  m <- newListArray ((1,1),(2,3)) [1,2,3,4,5,6] :: IO IOMatrix
  rowMeans m

array1dToUVectorD :: IO1dArray -> IO UVectorD
array1dToUVectorD array = do
  (=<<) (return . fromList) (getElems array)

smpdfs :: Int -> Int -> (UVectorD -> UVectorD) -> Int -> Int
       -> IO3dArray -> IO (Int, IO3dArray)
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
      dfmd = UV.foldr (+) 0 (UV.map abs fc)
  frthdf <- newArray ((1,1),(nd,nd+1)) 0 :: IO IOMatrix
  iejeitjtisjsls <- new 7 :: IO IOVectorI
  write iejeitjtisjsls 4 1
  write iejeitjtisjsls 5 2
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
                              let h = UV.map (*(2/(5*(fromInt nd +1))))
                                             (UV.zipWith (-) vi vj)
                                  ewd = UV.foldr (+) 0 (UV.map abs h)
                                  twoh = UV.map (*2) h
                                  t1 = f (UV.zipWith (-) cn twoh)
                                  t2 = f (UV.zipWith (+) cn twoh)
                                  t3 = UV.map (*6) fc
                                  t4 = f (UV.zipWith (-) cn h)
                                  t5 = f (UV.zipWith (+) cn h)
                                  t6 = UV.map (*(-4)) (UV.zipWith (+) t4 t5)
                                  tsum = (foldl1 (UV.zipWith (+))) [t1,t2,t3,t6]
                                  dfr1 = UV.foldr (+) 0 (UV.map abs tsum)
                                  dfr2 = if dfmd+dfr1/8 == dfmd then 0 else dfr1
                                  dfr3 = dfr2*ewd
                              dfmx <- UMV.read dfmxdfnx 0
                              case dfr3 >= dfmx of
                                True -> do
                                  is <- UMV.read iejeitjtisjsls 4
                                  js <- UMV.read iejeitjtisjsls 5
                                  write iejeitjtisjsls 2 is
                                  write iejeitjtisjsls 3 js
                                  write iejeitjtisjsls 4 i
                                  write iejeitjtisjsls 5 j
                                  write dfmxdfnx 1 dfmx
                                  write dfmxdfnx 0 dfr3
                                False -> do
                                  dfnx <- UMV.read dfmxdfnx 1
                                  case dfr3 >= dfnx of
                                    True -> do
                                      write iejeitjtisjsls 2 i
                                      write iejeitjtisjsls 3 j
                                      write dfmxdfnx 1 dfr3
                                    False -> return ()
                              writeArray frthdf (i,j) dfr3
                              case ewd >= y of
                                True -> do
                                  write iejeitjtisjsls 0 i
                                  write iejeitjtisjsls 1 j
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
          ie <- UMV.read iejeitjtisjsls 0
          je <- UMV.read iejeitjtisjsls 1
          write iejeitjtisjsls 4 ie
          write iejeitjtisjsls 5 je
        False -> do
          let loop :: Int -> Double -> Int -> IO Int
              loop l x ls | l == nd+2 = return ls
                          | otherwise = do
                            is <- UMV.read iejeitjtisjsls 4
                            js <- UMV.read iejeitjtisjsls 5
                            case (l /= is) && (l /= js) of
                              True -> do
                                let it = minimum [l,is,js]
                                    jt = maximum [l,is,js]
                                write iejeitjtisjsls 2 it
                                write iejeitjtisjsls 3 jt
                                let lt = is+js+l-it-jt
                                dfr1 <- readArray frthdf (it,lt)
                                dfr2 <- readArray frthdf (lt,jt)
                                let dfr = dfr1 + dfr2
                                case dfr >= x of
                                  True  -> loop (l+1) dfr l
                                  False -> loop (l+1) x ls
                              False -> loop (l+1) x ls
          ls <- loop 1 0 0
          write iejeitjtisjsls 6 ls
          is <- UMV.read iejeitjtisjsls 4
          js <- UMV.read iejeitjtisjsls 5
          difil <- readArray frthdf (min is ls, max is ls)
          diflj <- readArray frthdf (min js ls, max js ls)
          let dfnx = max difil diflj
          write dfmxdfnx 1 dfnx
          case dfmx/cuttb < dfnx && difil > diflj of
            True -> do
              is <- UMV.read iejeitjtisjsls 4
              js <- UMV.read iejeitjtisjsls 5
              it <- UMV.read iejeitjtisjsls 2
              write iejeitjtisjsls 2 is
              write iejeitjtisjsls 4 js
              write iejeitjtisjsls 5 it
            False -> return ()
  vrts2 <- newArray_ ((1,1,1),(nd,nd+1,sbs+nregions-1)) :: IO IO3dArray
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
  is <- UMV.read iejeitjtisjsls 4
  js <- UMV.read iejeitjtisjsls 5
  vtiIO <- extractColumn v is
  vtjIO <- extractColumn v js
  vti <- array1dToUVectorD vtiIO
  vtj <- array1dToUVectorD vtjIO
  case nregions == 4 of
    True -> do
      let vt = UV.map (/2) (UV.zipWith (+) vti vtj)
      replaceDimension vrts2 (js,top) vt
      replaceDimension vrts2 (is,sbs+1) vti
      replaceDimension vrts2 (js,sbs+1) vt
      replaceDimension vrts2 (is,sbs+2) vt
      replaceDimension vrts2 (is,sbs+3) vt
      replaceDimension vrts2 (js,sbs+2) vtj
      replaceDimension vrts2 (js,sbs+3) vtj
      it <- UMV.read iejeitjtisjsls 2
      jt <- UMV.read iejeitjtisjsls 3
      vtitIO <- mapIndices (1,nd) (\i -> (i,it,top)) vrts2 :: IO IO1dArray
      vtjtIO <- mapIndices (1,nd) (\i -> (i,jt,top)) vrts2 :: IO IO1dArray
      vtit <- (=<<) (return . fromList) (getElems vtitIO)
      vtjt <- (=<<) (return . fromList) (getElems vtjtIO)
      let vtt = UV.map (/2) (UV.zipWith (+) vtit vtjt)
      replaceDimension vrts2 (jt,top) vtt
      replaceDimension vrts2 (it,sbs+1) vtt
      replaceDimension vrts2 (jt,sbs+1) vtjt
      vti2IO <- mapIndices (1,nd) (\i -> (i,it,sbs+2)) vrts2 :: IO IO1dArray
      vtj2IO <- mapIndices (1,nd) (\i -> (i,jt,sbs+2)) vrts2 :: IO IO1dArray
      vti2 <- (=<<) (return . fromList) (getElems vti2IO)
      vtj2 <- (=<<) (return . fromList) (getElems vtj2IO)
      let vt2 = UV.map (/2) (UV.zipWith (+) vti2 vtj2)
      replaceDimension vrts2 (jt,sbs+2) vt2
      replaceDimension vrts2 (it,sbs+3) vt2
      replaceDimension vrts2 (jt,sbs+3) vtj2
    False -> do
      let vt = UV.map (/3) (UV.zipWith (+) (UV.map (*2) vti) vtj)
      replaceDimension vrts2 (js,top) vt
      replaceDimension vrts2 (is,sbs+1) vt
      case dfmx/cuttf < dfnx of
        True -> do
          replaceDimension vrts2 (js,sbs+1) vtj
          replaceDimension vrts2 (is,sbs+2) vt
          replaceDimension vrts2 (js,sbs+2) vtj
          ls <- UMV.read iejeitjtisjsls 6
          vtj1IO <- mapIndices (1,nd) (\i -> (i,js,sbs+1)) vrts2 :: IO IO1dArray
          vtl1IO <- mapIndices (1,nd) (\i -> (i,ls,sbs+1)) vrts2 :: IO IO1dArray
          vtj1 <- (=<<) (return . fromList) (getElems vtj1IO)
          vtl1 <- (=<<) (return . fromList) (getElems vtl1IO)
          let vt1 = UV.map (/2) (UV.zipWith (+) vtj1 vtl1)
          replaceDimension vrts2 (ls,sbs+1) vt1
          replaceDimension vrts2 (js,sbs+2) vt1
          replaceDimension vrts2 (ls,sbs+2) vtl1
        False -> do
          let vv = UV.map (/3) (UV.zipWith (+) vti (UV.map (*2) vtj))
          replaceDimension vrts2 (js,sbs+1) vv
          replaceDimension vrts2 (is,sbs+2) vv
          replaceDimension vrts2 (js,sbs+2) vtj
  return (nregions, vrts2)

replaceDimension :: IO3dArray -> (Int,Int) -> UVectorD -> IO ()
replaceDimension m (j,k) v = do
  (_, (n,_,_)) <- getBounds m
  let loop :: Int -> IO ()
      loop i | i == n+1 = return ()
             | otherwise = do
               writeArray m (i,j,k) ((UV.!) v (i-1))
               loop (i+1)
  loop 1

smpchc :: Int -> Int -> Int
smpchc nd key =
  if key == 0 || key == 3
    then div ((nd+4)*(nd+3)*(nd+2)) 6 + (nd+2)*(nd+1)
    else if key == 1
      then 2*nd + 3
      else if key == 2
        then div ((nd+3)*(nd+2)) 2 + 2*(nd+1)
        else div ((nd+5)*(nd+4)*(nd+3)*(nd+2)) 24 + 5*(div ((nd+2)*(nd+1)) 2)

check :: Int -> Int -> Int -> Double -> Double -> Int -> Int -> Int
check nd nf mxfs ea er sbs key =
  if ea < 0 && er < 0
    then 5
    else if nf < 1
      then 4
      else if nd < 2
        then 3
        else if key < 0 || key > 4
          then 2
          else if mxfs < sbs*smpchc nd key
            then 1
            else 0

adsimp :: Int -> Int -> Int -> (UVectorD -> UVectorD) -> Double -> Double
       -> Int -> IO3dArray -> Bool -> IO (UVectorD, UVectorD, Int, Bool)
adsimp nd nf mxfs f ea er key vrts info = do
  case key == 0 of
    True -> adsimp nd nf mxfs f ea er 3 vrts info
    False -> do
      (_, (_,_,sbs)) <- getBounds vrts
      let b = smpchc nd key
      smpsad nd nf f mxfs ea er key b sbs vrts info

smpsad_old :: Int -> Int -> (UVectorD -> UVectorD) -> Int -> Double -> Double -> Int
       -> Int -> Int -> IO3dArray -> Bool -> IO (UVectorD, UVectorD, Int, Bool)
smpsad_old nd nf f mxfs ea er key rcls sbs vrts info = do
  let dfcost = 1 + 2*nd*(nd+1)
  -- vl <- newArray (1,nf) 0 :: IO IO1dArray
  -- ae <- newArray (1,nf) 0 :: IO IO1dArray
  (g, w, ptsIO) <- smprms nd key
  pts <- freeze ptsIO
  -- let fct = fromInt (product [1..nd]) -- je m'en sers ?
  -- vls <- newArray_ ((1,1),(nf,sbs)) :: IO IOMatrix
  -- aes <- newArray_ ((1,1),(nf,sbs)) :: IO IOMatrix
  simplices <- mapM (toSimplex vrts (nd+1)) [1..sbs]
  let vol = S.fromList $ map simplexVolume simplices
      nv = sbs*rcls
  matrices <- mapM
              (\k -> mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,k)) vrts)
              [1..sbs]
  (basvals, rgnerrs) <- (=<<) (return.unzip)
                        (mapM (\k ->
                                smprul (matrices!!k) nf f (index vol k) g w pts)
                        [0..(sbs-1)])
  -- aes <- matrixFromColumns rgnerrs
  -- vls <- matrixFromColumns basvals
  rgnerrs' <- mapM array1dToUVectorD rgnerrs
  basvals' <- mapM array1dToUVectorD basvals
  let vl = foldl1 (UV.zipWith (+)) basvals'
      ae = foldl1 (UV.zipWith (+)) rgnerrs'
      aes = V.fromListN sbs rgnerrs'
      vls = V.fromListN sbs basvals'
  let fl = UV.any (> (max ea (UV.maximum (UV.map ((*er).abs) vl)))) ae
  -- let loop :: Params -> IO (UVectorD, UVectorD, Int, Bool)
  --     loop params | not fl || (nv+dfcost+4*rcls > mxfs) = return (vl, ae, nv, fl)
  let loop :: Params -> IO (UVectorD, UVectorD, Int, Bool)
      loop params | not fl || (nv+dfcost+4*rcls > mxfs) = return (vl, ae, nv, fl) --(rowSumsV vls, rowSumsV aes, nv, fl)
                  | otherwise = do
                    let maxs = V.map (UV.maximum) aes
                        id = fromJust $ V.findIndex (== V.maximum maxs) maxs
                        vl0 = UV.zipWith (-) vl (vls V.! id)
                        ae0 = UV.zipWith (-) ae (aes V.! id)
                    (nregions, vrts2) <- smpdfs nd nf f (id+1) sbs vrts
                    let vi = (index vol id)/(fromInt nregions)
                        nv2 = nv + (nregions-1)*rcls -- nregions*rcls ?
                        sbs2 = sbs + (nregions-1)
                    matrices <- mapM
                                (\k -> mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,k)) vrts2)
                                ((id+1):[(sbs+1)..(sbs+nregions-1)])
                    (basvals, rgnerrs) <- (=<<) (return.unzip)
                                          (mapM (\m ->
                                                  smprul m nf f vi g w pts)
                                          matrices)
                    -- TODO améliorer cette fonction - on a besoin que des rowSums
                    -- changer vol dans smprul n'est qu'un facteur je pense
                    -- je pense que vl et ae sont ok, et pas besoin de aes vls
                    -- sauf aes pour les max => garder juste les max ?
                    -- oui c'est bon !
                    -- reste just à améliorer vol2
                    -- est-ce qu'un vol mutable serait le plus efficace ?
                    -- peut-être bof pour faire ++ replicate
                    rgnerrs' <- mapM array1dToUVectorD rgnerrs
                    basvals' <- mapM array1dToUVectorD basvals
                    let vl2 = UV.zipWith (+) vl0 (foldl1 (UV.zipWith (+)) basvals')
                        ae2 = UV.zipWith (+) ae0 (foldl1 (UV.zipWith (+)) rgnerrs')
                        aes2 = (aes V.// [(id, head rgnerrs')]) V.++ (V.tail (V.fromList rgnerrs')) -- fromListN (nregions-1)
                        vls2 = (vls V.// [(id, head basvals')]) V.++ (V.tail (V.fromList basvals'))
                        fl2 = UV.any (> (max ea (UV.maximum (UV.map ((*er).abs) vl2)))) ae2
                        vol2 = (update id vi vol) >< (S.replicate (nregions-1) vi)
                        -- vol2_temp = vol ++ (replicate (nregions-1) vi)
                        -- vol2 = map (\i -> if (i == id) then vi else (vol2_temp!!i)) [0..(length vol2_temp - 1)]
                    loop (fl2, nv2, sbs2, aes2, vls2, ae2, vl2, vrts2, vol2)
                  where
                    (fl, nv, sbs, aes, vls, ae, vl, vrts, vol) = params
  -- dans le code R il fait rowSums mais ça me semble inutile
  loop (fl, nv, sbs, aes, vls, ae, vl, vrts, vol)

type Params =
  (Bool, Int, Int, V.Vector UVectorD, V.Vector UVectorD, UVectorD, UVectorD, IO3dArray, Seq Double)
--
type Params2 =
  (Bool, Int, Int, Seq UVectorD, Seq UVectorD, UVectorD, UVectorD, IO3dArray, Seq Double)

smpsad :: Int -> Int -> (UVectorD -> UVectorD) -> Int -> Double -> Double -> Int
       -> Int -> Int -> IO3dArray -> Bool -> IO (UVectorD, UVectorD, Int, Bool)
smpsad nd nf f mxfs ea er key rcls sbs vrts info = do
  let dfcost = 1 + 2*nd*(nd+1)
  (g, w, ptsIO) <- smprms nd key
  pts <- freeze ptsIO
  simplices <- mapM (toSimplex vrts (nd+1)) [1..sbs]
  let vol = S.fromList $ map simplexVolume simplices
      nv = sbs*rcls
  matrices <- mapM
              (\k -> mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,k)) vrts)
              (S.fromList [1..sbs])
  br <- (mapM (\k -> smprul (index matrices k) nf f (index vol k) g w pts)
                        (S.fromList [0..(sbs-1)]))
  rgnerrs' <- mapM (array1dToUVectorD.snd) br
  basvals' <- mapM (array1dToUVectorD.fst) br
  let vl = foldl1 (UV.zipWith (+)) basvals'
      ae = foldl1 (UV.zipWith (+)) rgnerrs'
      aes = rgnerrs' -- V.fromList (DF.toList rgnerrs')
      vls = basvals' -- V.fromList (DF.toList basvals')
  let fl = UV.any (> (max ea (UV.maximum (UV.map ((*er).abs) vl)))) ae
  let loop :: Params2 -> IO (UVectorD, UVectorD, Int, Bool)
      loop params | not fl || (nv+dfcost+4*rcls > mxfs) = return (vl, ae, nv, fl) --(rowSumsV vls, rowSumsV aes, nv, fl)
                  | otherwise = do
                    let maxs = fmap UV.maximum aes
                        id = fromJust $ S.findIndexL (== maximum maxs) maxs
                        vl0 = UV.zipWith (-) vl (index vls id)
                        ae0 = UV.zipWith (-) ae (index aes id)
                    (nregions, vrts2) <- smpdfs nd nf f (id+1) sbs vrts
                    let vi = (index vol id)/(fromInt nregions)
                        nv2 = nv + (nregions-1)*rcls -- nregions*rcls ?
                        sbs2 = sbs + (nregions-1)
                    matrices <- mapM
                                (\k -> mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,k)) vrts2)
                                (S.fromList ((id+1):[(sbs+1)..(sbs+nregions-1)]))
                    br <- (mapM (\m -> smprul m nf f vi g w pts)
                                          matrices)
                    rgnerrs' <- mapM (array1dToUVectorD.snd) br
                    basvals' <- mapM (array1dToUVectorD.fst) br
                    let vl2 = UV.zipWith (+) vl0 (foldl1 (UV.zipWith (+)) basvals')
                        ae2 = UV.zipWith (+) ae0 (foldl1 (UV.zipWith (+)) rgnerrs')
                        aes2 = (update id (index rgnerrs' 0) aes) >< (S.drop 1 rgnerrs')
                        vls2 = (update id (index basvals' 0) vls) >< (S.drop 1 basvals')
                        fl2 = UV.any (> (max ea (UV.maximum (UV.map ((*er).abs) vl2)))) ae2
                        vol2 = (update id vi vol) >< (S.replicate (nregions-1) vi)
                        -- vol2_temp = vol ++ (replicate (nregions-1) vi)
                        -- vol2 = map (\i -> if (i == id) then vi else (vol2_temp!!i)) [0..(length vol2_temp - 1)]
                    loop (fl2, nv2, sbs2, aes2, vls2, ae2, vl2, vrts2, vol2)
                  where
                    (fl, nv, sbs, aes, vls, ae, vl, vrts, vol) = params
  -- dans le code R il fait rowSums mais ça me semble inutile
  loop (fl, nv, sbs, aes, vls, ae, vl, vrts, vol)

rowSumsV :: V.Vector UVectorD -> UVectorD
rowSumsV = V.foldl1 (UV.zipWith (+))

matrixFromColumns :: [IO1dArray] -> IO UMatrix
matrixFromColumns columns = do
  (_,nrow) <- getBounds (head columns)
  let ncol = length columns
  u1darrays <- mapM IOA.freeze columns :: IO [UArray Int Double]
  let assocs = map UA.assocs u1darrays
  let assocs' = map (\k -> map (\(i,e) -> ((k,i),e)) (assocs !! (k-1))) [1..ncol]
  return $ UA.array ((1,1),(nrow,ncol)) (concat assocs')


toSimplex :: IO3dArray -> Int -> Int -> IO Simplex
toSimplex m n k = do
  (_, (nrow,_,_)) <- getBounds m
  let getColumn :: Int -> IO IO1dArray
      getColumn col = mapIndices (1,nrow) (\i -> (i,col,k)) m
  columns <- mapM getColumn [1..n]
  mapM getElems columns
