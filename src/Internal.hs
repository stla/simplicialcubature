{-# LANGUAGE BangPatterns #-}
module Internal
  where
import           Common
import           Data.Array.IArray           (Array)
import           Data.Array.IO               (IOUArray, getBounds, getElems,
                                              mapIndices, newArray, newArray_,
                                              newListArray, readArray,
                                              writeArray)
import           Data.List                   (foldl', foldl1')
import           Control.Monad               ((<$!>), (=<<))
import           Data.Array.Unboxed          (UArray)
import qualified Data.Foldable               as DF
import           Data.Maybe                  (fromJust)
import           Data.Sequence               (Seq, index, update, (><))
import qualified Data.Sequence               as S
import           Data.Vector.Unboxed         (Vector, fromList, toList,
                                              unsafeFreeze)
import qualified Data.Vector.Unboxed         as UV
import           Data.Vector.Unboxed.Mutable (IOVector, new, unsafeRead,
                                              unsafeWrite)
import qualified Data.Vector.Unboxed.Mutable as UMV
import           Math.Combinat.Permutations  (permuteMultiset)
import           Simplex                     (Simplex, simplexVolume)

type IOMatrix = IOUArray (Int,Int) Double
type IO1dArray = IOUArray Int Double
type IMatrix = Array (Int,Int) Double
type UMatrix = UArray (Int,Int) Double
type IOVectorD = IOVector Double
type IOVectorI = IOVector Int
type UVectorD = Vector Double
type UVectorI = Vector Int

toDbl :: Int -> Double
toDbl = fromIntegral

pow :: Double -> Int -> Double
pow x n = foldr (*) 1 (replicate n x)

square :: Double -> Double
square x = x*x

smprms :: Int -> Int -> IO (IOMatrix, Seq UVectorD, [Int])
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
      ndbl = toDbl n
      sqrt15 = 3.872983346207416885179265399782399611
      r1 = (ndbl + 4 - sqrt15) / (ndbl*ndbl + 8*ndbl + 1)
      s1 = 1 - ndbl*r1
      l1 = s1 - r1
  _ <- mapM (\j -> writeArray g (j,1) (1/(ndbl+1))) [1 .. np]
  unsafeWrite pts 0 1
  writeArray g (1,gms+1) s1
  _ <- mapM (\j -> writeArray g (j,gms+1) r1) [2 .. np]
  unsafeWrite pts gms np
  case key < 4 of
    True -> do
      writeArray w (1,rls) 1
      writeArray w (gms+1,rls-1) (1/(ndbl+1))
    False -> return ()
  let iw = if key < 4 then rls-2 else rls
  writeArray g (1,2) (3/(ndbl+3))
  _ <- mapM (\j -> writeArray g (j,2) (1/(ndbl+3))) [2 .. np]
  unsafeWrite pts 1 np
  let n2double = toDbl n2
  writeArray w (2,iw) (pow (ndbl+3) 2 / (4*n2double))
  case key > 1 of
    True -> do
      case n == 2 of
        True -> do
          let l2 = 0.62054648267200632589046034361711
              r1' = (1+sqrt(0.5 - l2*l2))/3
          writeArray g (1,gms+1) (1 - 2*r1')
          _ <- mapM (\j -> writeArray g (j,gms+1) r1') [2 .. np]
          unsafeWrite pts gms 3
          writeArray w (gms+1,iw-1) (1/6)
          let r2 = (1-l2)/3
              s2 = 1 - 2*r2
          writeArray g (1,gms+2) s2
          _ <- mapM (\j -> writeArray g (j,gms+2) r2) [2 .. np]
          unsafeWrite pts (gms+1) 3
          writeArray w (gms+2,iw-1) (1/6)
        False -> do
          let r2 = (ndbl+4+sqrt15) / (ndbl*ndbl+8*ndbl+1)
              s2 = 1 - ndbl*r2
              l2 = s2 - r2
          writeArray g (1,gms+2) s2
          _ <- mapM (\j -> writeArray g (j,gms+2) r2) [2 .. np]
          unsafeWrite pts (gms+1) np
          writeArray w (gms+2,iw-1) ((2/(ndbl+3)-l1)/(n2double*(l2-l1)*l2*l2))
          writeArray w (gms+1,iw-1) ((2/(ndbl+3)-l2)/(n2double*(l1-l2)*l1*l1))
      writeArray g (1,3) (5/(ndbl+5))
      _ <- mapM (\j -> writeArray g (j,3) (1/(ndbl+5))) [2 .. np]
      unsafeWrite pts 2 np
      writeArray g (1,4) (3/(ndbl+5))
      writeArray g (2,4) (3/(ndbl+5))
      _ <- mapM (\j -> writeArray g (j,4) (1/(ndbl+5))) [3 .. np]
      unsafeWrite pts 3 (div (n*np) 2)
      let n4double = toDbl n4
      writeArray w (2,iw-2) (-(pow (ndbl+3) 5) / (16*n4double))
      writeArray w (3,iw-2) (pow (ndbl+5) 5 / (16*n4double*(ndbl+5)))
      writeArray w (4,iw-2) (pow (ndbl+5) 5 / (16*n4double*(ndbl+5)))
    False -> return ()
  case key > 2 of
    True -> do
      let u1 = (ndbl+7+2*sqrt15) / (ndbl*ndbl+14*ndbl-11)
          v1 = (1-(ndbl-1)*u1)/2
          d1 = v1 - u1
      writeArray g (1,gms+3) v1
      writeArray g (2,gms+3) v1
      _ <- mapM (\j -> writeArray g (j,gms+3) u1) [3 .. np]
      unsafeWrite pts (gms+2) (div (n*np) 2)
      let u2 = (ndbl+7-2*sqrt15) / (ndbl*ndbl+14*ndbl-11)
          v2 = (1-(ndbl-1)*u2)/2
          d2 = v2 - u2
      writeArray g (1,gms+4) v2
      writeArray g (2,gms+4) v2
      _ <- mapM (\j -> writeArray g (j,gms+4) u2) [3 .. np]
      unsafeWrite pts (gms+3) (div (n*np) 2)
      case n == 2 of
        True -> do
          writeArray w (gms+3,iw-3) ((155-sqrt15)/1200)
          writeArray w (gms+4,iw-3) ((155+sqrt15)/1200)
          writeArray w (1,iw-3) (1 - 31/40)
        False -> do
          case n == 3 of
            True -> do
              writeArray w (gms+1,iw-3) ((2665+14*sqrt15)/37800)
              writeArray w (gms+2,iw-3) ((2665-14*sqrt15)/37800)
              writeArray w (gms+3,iw-3) (30/567)
              unsafeWrite pts (gms+3) 0
            False -> do
              let r2 = (ndbl+4+sqrt15) / (ndbl*ndbl+8*ndbl+1)
                  l2 = 1 - (ndbl+1)*r2
                  den = pow l1 4 * (l1-l2) * (toDbl n4)
              writeArray w (gms+1,iw-3)
                           ((2*(27-ndbl)/(ndbl+5)-l2*(13-ndbl))/den)
              writeArray w (gms+2,iw-3)
                           ((2*(27-ndbl)/(ndbl+5)-l1*(13-ndbl))/den)
              writeArray w (gms+3,iw-3)
                           ((2/(ndbl+5)-d2)/((toDbl n4)*(d1-d2)*(pow d1 4)))
              writeArray w (gms+4,iw-3)
                           ((2/(ndbl+5)-d1)/((toDbl n4)*(d2-d1)*(pow d2 4)))
      writeArray g (1,5) (7/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,5) (1/(ndbl+7))) [2 .. np]
      unsafeWrite pts 4 np
      writeArray g (1,6) (5/(ndbl+7))
      writeArray g (2,6) (3/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,6) (1/(ndbl+7))) [3 .. np]
      unsafeWrite pts 5 (np*n)
      _ <- mapM (\i -> writeArray g (i,7) (3/(ndbl+7))) [1,2,3]
      _ <- case np > 3 of
        True  -> mapM (\i -> writeArray g (i,7) (1/(ndbl+7))) [4..np]
        False -> return [()]
      unsafeWrite pts 6 (div ((n-1)*n*np) 6)
      writeArray w (2,iw-4) ((pow (ndbl+3) 7)/(toDbl $ 128*n4*(n+5)))
      _ <- mapM
           (\i -> writeArray w (i,iw-4) (-(pow (ndbl+5) 7) /(toDbl $ 64*n6)))
           [3,4]
      _ <- mapM
           (\i -> writeArray w (i,iw-4) (pow (ndbl+7) 7/(toDbl $ 64*n6*(n+7))))
           [5,6,7]
      return ()
    False -> return ()
  case key == 4 of
    True -> do
      let sg = 1/(toDbl $ 23328*n6)
          u5 = -216 * sg * (toDbl $ 52212 - n*(6353 + n*(1934-n*27)))
          u6 =  1296 * sg * (toDbl $ 7884 - n*(1541 - n*9))
          u7 = -7776 * sg * (toDbl $ 8292 - n*(1139 - n*3))/(ndbl + 7)
          p0 = -144 * (142528 + n*(23073 - n*115))
          p1 = -12 * (6690556 + n*(2641189 + n*(245378 - n*1495)))
          p2 = -16 * (6503401 + n*(4020794+n*(787281+n*(47323-n*385))))
          p3 = -(6386660 + n*(4411997+n*(951821+n*(61659-n*665))))*(n + 7)
          a = (toDbl p2)/(toDbl $ 3*p3)
          p = a*((toDbl p1)/(toDbl p2) - a)
          q = a*(2*a*a - (toDbl p1)/(toDbl p3)) + (toDbl p0)/(toDbl p3)
          th = acos(-q/(2*(sqrt(-p*p*p))))/3
          r = 2*sqrt(-p)
          tp = 2*pi/3
          a1 = -a + r*(cos(th))
          a2 = -a + r*(cos(th+2*tp))
          a3 = -a + r*(cos(th+tp))
          npdbl = toDbl np
      writeArray g (1,gms+5) ((1-ndbl*a1)/npdbl)
      _ <- mapM (\i -> writeArray g (i,gms+5) ((1+a1)/npdbl)) [2..np]
      unsafeWrite pts (gms+4) np
      writeArray g (1,gms+6) ((1-ndbl*a2)/npdbl)
      _ <- mapM (\i -> writeArray g (i,gms+6) ((1+a2)/npdbl)) [2..np]
      unsafeWrite pts (gms+5) np
      writeArray g (1,gms+7) ((1-ndbl*a3)/npdbl)
      _ <- mapM (\i -> writeArray g (i,gms+7) ((1+a3)/npdbl)) [2..np]
      unsafeWrite pts (gms+6) np
      writeArray w (gms+5,iw-5)
                 ((u7-(a2+a3)*u6+a2*a3*u5)/(a1*a1-(a2+a3)*a1+a2*a3)/(pow a1 5))
      writeArray w (gms+6,iw-5)
                 ((u7-(a1+a3)*u6+a1*a3*u5)/(a2*a2-(a1+a3)*a2+a1*a3)/(pow a3 5))
      writeArray w (gms+7,iw-5)
                 ((u7-(a2+a1)*u6+a2*a1*u5)/(a3*a3-(a2+a1)*a3+a2*a1)/(pow a3 5))
      writeArray g (1,gms+8) (4/(ndbl+7))
      writeArray g (2,gms+8) (4/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,gms+8) (1/(ndbl+7))) [3..np]
      unsafeWrite pts (gms+7) (div (np*n) 2)
      writeArray w (gms+8,iw-5) (10*(pow (ndbl+7) 6)/(toDbl $ 729*n6))
      writeArray g (1,gms+9) (11/(ndbl+7)/2)
      writeArray g (2,gms+9) (5/(ndbl+7)/2)
      _ <- mapM (\i -> writeArray g (i,gms+9) (1/(ndbl+7))) [3..np]
      unsafeWrite pts (gms+8) (np*n)
      writeArray w (gms+9,iw-5) (64*(pow (ndbl+7) 6) / (toDbl $ 6561*n6))
      writeArray w (4,iw-5) (-(pow (ndbl+5) 7) / (toDbl $ 64*n6))
      writeArray w (7,iw-5) (pow (ndbl+7) 7 / (toDbl $ 64*n6*(n+7)))
      writeArray g (1,8) (9/(ndbl+9))
      _ <- mapM (\i -> writeArray g (i,8) (1/(ndbl+9))) [2..np]
      unsafeWrite pts 7 np
      writeArray g (1,9) (7/(ndbl+9))
      writeArray g (2,9) (3/(ndbl+9))
      _ <- mapM (\i -> writeArray g (i,9) (1/(ndbl+9))) [3..np]
      unsafeWrite pts 8 (np*n)
      _ <- mapM (\i -> writeArray g (i,10) (5/(ndbl+9))) [1,2]
      _ <- mapM (\i -> writeArray g (i,10) (1/(ndbl+9))) [3..np]
      unsafeWrite pts 9 (div (np*n) 2)
      writeArray g (1,11) (5/(ndbl+9))
      _ <- mapM (\i -> writeArray g (i,11) (3/(ndbl+9))) [2,3]
      _ <- case np > 3 of
        True  -> mapM (\i -> writeArray g (i,11) (1/(ndbl+9))) [4..np]
        False -> return [()]
      unsafeWrite pts 10 (div (np*n*(n-1)) 2)
      writeArray w (2,iw-6) (-(pow (ndbl+3) 9)/(toDbl $ 6*256*n6))
      _ <- mapM (\i -> writeArray w (i,iw-6)
                                  (pow (ndbl+5) 9 /(toDbl $ 512*n6*(n+7))))
                [3,4]
      _ <- mapM (\i -> writeArray w (i,iw-6)
                                  (-(pow (ndbl+7) 9)/(toDbl $ 256*n8)))
                [5,6,7]
      _ <- mapM (\i -> writeArray w (i,iw-6)
                                  (pow (ndbl+9) 9/(toDbl $ 256*n8*(n+9))))
                [8..11]
      case n > 2 of
        True ->  do
          _ <- mapM (\i -> writeArray g (i,12) (3/(ndbl+9))) [1..4]
          _ <- mapM (\i -> writeArray g (i,12) (1/(ndbl+9))) [5..np]
          unsafeWrite pts 11 (div (np*n*(n-1)*(n-2)) 24)
          writeArray w (12,iw-6) (pow (ndbl+9) 9 / (toDbl $ 256*n8*(n+9)))
        False -> return ()
    False -> return ()
  rowsIO <- mapM (extractRow w) (S.fromList [2..wts])
  rows <- mapM array1dToUVectorD rowsIO
  let cols = transpose rows
  pts_out <- UV.unsafeFreeze pts
  let ptsU = UV.map toDbl pts_out
      row1 = fmap (\col -> 1 - UV.foldr (+) 0
                           (UV.zipWith (*) (UV.tail ptsU) col))
                  (S.take rls cols)
      wcols = fmap (\j -> UV.cons (index row1 j) (index cols j))
                   (S.fromList [0..(rls-1)])
      col1 = index wcols 0
      wcols2 = (S.<|) col1
                      (fmap (\col -> UV.zipWith (-) col col1) (S.drop 1 wcols))
      col2 = index wcols2 1
      nb = UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map square col1))
      ratio = nb / (UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map square col2)))
      wcol2 = UV.map (*(sqrt ratio)) col2
      wcols3 = (S.<|) col1 ((S.<|) wcol2 (S.drop 2 wcols2))
  let updateW :: Seq UVectorD -> Int -> Seq UVectorD
      updateW cols' k = update (k-1) wknew cols'
         where
          ptsW = UV.map (/nb) (UV.zipWith (*) ptsU (index cols' (k-1)))
          slice = S.drop 1 (S.take (k-1) cols')
          prod1 = (fromList.(DF.toList)) $ -- ou alors pas de vecteurs
                  fmap ((UV.foldr (+) 0).(UV.zipWith (*) ptsW)) slice
          rows' = transpose slice
          prod2 = (fromList.(DF.toList)) $
                  fmap ((UV.foldr (+) 0).(UV.zipWith (*) prod1)) rows'
          wk = UV.zipWith (-) (index cols' (k-1)) prod2
          ratio' = nb / (UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map square wk)))
          wknew = UV.map (*(sqrt ratio')) wk
      wcolsnew = foldl' updateW wcols3 [3..rls]
  return (g, transpose wcolsnew, toList $ UV.findIndices (/= 0) pts_out)

transpose :: Seq UVectorD -> Seq UVectorD
transpose cols =
  fmap (\i -> (fromList.(DF.toList)) $ fmap (\col -> col UV.!i) cols)
       (S.fromList ([0..(UV.length (index cols 0) - 1)]))

matprod :: IOMatrix -> UVectorD -> IO UVectorD
matprod mat x = do
  (_, (m,n)) <- getBounds mat
  out <- UMV.new m :: IO IOVectorD
  let step i | i == m+1 = unsafeFreeze out
             | otherwise = do
              !coef <- innerstep 1 0
              unsafeWrite out (i-1) coef
              step (i+1)
              where
                innerstep :: Int -> Double -> IO Double
                innerstep j !s | j == n+1 = return s
                               | otherwise = do
                                 mat_ij <- readArray mat (i,j)
                                 innerstep (j+1) (s + mat_ij * (x UV.! (j-1)))
  step 1

smpsms :: IOMatrix -> Int -> (UVectorD -> UVectorD) -> IO1dArray
       -> Double -> IO IO1dArray
smpsms vertex nf f g scalar = do
  gAsList <- getElems g
  f_gPermuts <- mapM ((fmap f).(matprod vertex).fromList)
                     (permuteMultiset gAsList)
  newListArray (1,nf)
               (toList (UV.map (*scalar) (foldl1' (UV.zipWith (+)) f_gPermuts)))

extractColumn :: IOMatrix -> Int -> IO IO1dArray
extractColumn m j = do
  (_, (nrow,_)) <- getBounds m
  mapIndices (1,nrow) (\i -> (i,j)) m

extractRow :: IOMatrix -> Int -> IO IO1dArray
extractRow m i = do
  (_, (_,ncol)) <- getBounds m
  mapIndices (1,ncol) (\j -> (i,j)) m

outerProduct :: IO1dArray -> UVectorD -> IO IOMatrix
outerProduct x1 x2 = do
  (_, n1) <- getBounds x1
  let n2 = UV.length x2
  out <- newArray_ ((1,1),(n1,n2)) :: IO IOMatrix
  let step :: Int -> IO IOMatrix
      step i | i == n1+1 = return out
             | otherwise = do
                x1_i <- readArray x1 i
                inner x1_i 0
              where
                inner !x j | j == n2 = step (i+1)
                           | otherwise = do
                              writeArray out (i,j+1) (x * (x2 UV.! j))
                              inner x (j+1)
  step 1

sumMatrices :: [IOMatrix] -> IO IOMatrix
sumMatrices matrices = do
  (_, (n1,n2)) <- getBounds (head matrices)
  out <- newArray_ ((1,1),(n1,n2)) :: IO IOMatrix
  let step :: Int -> IO IOMatrix
      step i | i == n1+1 = return out
             | otherwise = inner 1
               where
                 inner :: Int -> IO IOMatrix
                 inner j | j == n2+1 = step (i+1)
                         | otherwise = do
                           coefs <- mapM (\m -> readArray m (i,j)) matrices
                           writeArray out (i,j) (sum coefs)
                           inner (j+1)
  step 1

smprul :: IOMatrix -> Int -> (UVectorD -> UVectorD) -> Double -> IOMatrix
       -> Seq UVectorD -> [Int] -> IO (IO1dArray, IO1dArray)
smprul vrts nf f vol g w pospts = do
  let rtmn = 0.1
      small = 1e-12
      errcof = 8
      rls = UV.length (index w 0)
--      ptsPositive = toList $ UV.findIndices (> 0) pts
  toSum <- mapM (\k -> do
                         g_colk <- extractColumn g (k+1)
                         sms <- smpsms vrts nf f g_colk vol
                         outerProduct sms (index w k))
                pospts
  rule <- sumMatrices toSum
  basval <- extractColumn rule 1
  rgnerr <- newArray (1,nf) 0 :: IO IO1dArray
  let step :: Int -> IO ()
      step i | i == nf+1 = return ()
             | otherwise = do
                basval_i <- readArray basval i
                let nmbs = abs basval_i
                !(rt, nmcp) <- innerstep rls rtmn nmbs
                case rt < 1 && rls > 3 of
                  True  -> writeArray rgnerr i (rt*nmcp)
                  False -> return ()
                rgnerr_i <- readArray rgnerr i
                writeArray rgnerr i (max (errcof*rgnerr_i) (small*nmbs))
                step (i+1)
              where
               innerstep :: Int -> Double -> Double -> IO (Double, Double)
               innerstep k !x !y | k == 1 = return (x, y)
                                 | otherwise = do
                                  rule_ik <- readArray rule (i,k)
                                  rule_ikm1 <- readArray rule (i,k-1)
                                  let nmrl = max (abs rule_ik) (abs rule_ikm1)
                                  rgnerr_i <- readArray rgnerr i
                                  writeArray rgnerr i (max nmrl rgnerr_i)
                                  case nmrl > small*y && k < rls of
                                    True -> innerstep (k-2) (max (nmrl/y) x)
                                                      nmrl
                                    False -> innerstep (k-2) x nmrl
  step 1
  return (basval, rgnerr)

rowMeans :: IOMatrix -> IO UVectorD
rowMeans m = do
  (_, (nrow,ncol)) <- getBounds m
  outIO <- new nrow :: IO IOVectorD
  let step :: Int -> IO ()
      step i | i == nrow+1 = return ()
             | otherwise = do
                !sum_i <- inner 1 0
                unsafeWrite outIO (i-1) sum_i
                step (i+1)
              where
                inner :: Int -> Double -> IO Double
                inner j !x | j == ncol+1 = return (x / (toDbl ncol))
                           | otherwise = do
                             coef <- readArray m (i,j)
                             inner (j+1) (x + coef)
  step 1
  unsafeFreeze outIO

array1dToUVectorD :: IO1dArray -> IO UVectorD
array1dToUVectorD array = (<$!>) fromList (getElems array)

getVectors :: Int -> IO3dArray -> Int -> Int -> Int -> IO (UVectorD, UVectorD) -- pas de gain
getVectors n m k j1 j2 = do
  out1 <- new n :: IO IOVectorD
  out2 <- new n :: IO IOVectorD
  let loop :: Int -> IO ()
      loop i | i == n+1 = return ()
             | otherwise = do
               coef1 <- readArray m (i,j1,k)
               unsafeWrite out1 (i-1) coef1
               coef2 <- readArray m (i,j2,k)
               unsafeWrite out2 (i-1) coef2
               loop (i+1)
  loop 1
  out1U <- UV.unsafeFreeze out1
  out2U <- UV.unsafeFreeze out2
  return (out1U,out2U)


smpdfs :: Int -> (UVectorD -> UVectorD) -> Int -> Int
       -> IO3dArray -> IO (Int, IO3dArray)
smpdfs nd f top sbs vrts = do
  let cuttf = 2.0
      cuttb = 8.0
  v <- mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,top)) vrts
  cn <- rowMeans v
  let fc = f cn
      dfmd = UV.foldr (+) 0 (UV.map abs fc)
  frthdf <- newArray ((1,1),(nd,nd+1)) 0 :: IO IOMatrix
  iejeitjtisjsls <- new 7 :: IO IOVectorI
  unsafeWrite iejeitjtisjsls 4 1
  unsafeWrite iejeitjtisjsls 5 2
  dfmxdfnx <- UMV.replicate 2 0 :: IO IOVectorD
  let step :: Int -> Double -> IO ()
      step i x | i == nd+1 = return ()
               | otherwise = do
                  !emx <- inner (i+1) x
                  step (i+1) emx
                where
                  inner :: Int -> Double -> IO Double
                  inner j !y | j == nd+2 = return y
                             | otherwise = do
                              vi <- (=<<) array1dToUVectorD (extractColumn v i)
                              vj <- (=<<) array1dToUVectorD (extractColumn v j)
                              let h = UV.map (*(2/(5*(toDbl nd +1))))
                                             (UV.zipWith (-) vi vj)
                                  ewd = UV.foldr (+) 0 (UV.map abs h)
                                  twoh = UV.map (*2) h
                                  t1 = f (UV.zipWith (-) cn twoh)
                                  t2 = f (UV.zipWith (+) cn twoh)
                                  t3 = UV.map (*6) fc
                                  t4 = f (UV.zipWith (-) cn h)
                                  t5 = f (UV.zipWith (+) cn h)
                                  t6 = UV.map (*(-4)) (UV.zipWith (+) t4 t5)
                                  tsum = foldl1' (UV.zipWith (+)) [t1,t2,t3,t6]
                                  dfr1 = UV.foldr (+) 0 (UV.map abs tsum)
                                  dfr2 = if dfmd+dfr1/8 == dfmd then 0 else dfr1
                                  dfr3 = dfr2*ewd
                              dfmx <- unsafeRead dfmxdfnx 0
                              case dfr3 >= dfmx of
                                True -> do
                                 is <- unsafeRead iejeitjtisjsls 4
                                 js <- unsafeRead iejeitjtisjsls 5
                                 unsafeWrite iejeitjtisjsls 2 is
                                 unsafeWrite iejeitjtisjsls 3 js
                                 unsafeWrite iejeitjtisjsls 4 i
                                 unsafeWrite iejeitjtisjsls 5 j
                                 unsafeWrite dfmxdfnx 1 dfmx
                                 unsafeWrite dfmxdfnx 0 dfr3
                                False -> do
                                 dfnx <- unsafeRead dfmxdfnx 1
                                 case dfr3 >= dfnx of
                                  True -> do
                                    unsafeWrite iejeitjtisjsls 2 i
                                    unsafeWrite iejeitjtisjsls 3 j
                                    unsafeWrite dfmxdfnx 1 dfr3
                                  False -> return ()
                              writeArray frthdf (i,j) dfr3
                              case ewd >= y of
                                True -> do
                                  unsafeWrite iejeitjtisjsls 0 i
                                  unsafeWrite iejeitjtisjsls 1 j
                                  inner (j+1) ewd
                                False -> inner (j+1) y
  step 1 0
  dfmx <- unsafeRead dfmxdfnx 0
  dfnx <- unsafeRead dfmxdfnx 1
  let nregions = if dfnx > dfmx/cuttf then 4 else 3
  case dfnx > dfmx/cuttf of
    True -> return ()
    False -> do
      case dfmx == 0 of
        True -> do
          ie <- unsafeRead iejeitjtisjsls 0
          je <- unsafeRead iejeitjtisjsls 1
          unsafeWrite iejeitjtisjsls 4 ie
          unsafeWrite iejeitjtisjsls 5 je
        False -> do
          let loop :: Int -> Double -> Int -> IO Int
              loop l !x !ls | l == nd+2 = return ls
                            | otherwise = do
                              is <- unsafeRead iejeitjtisjsls 4
                              js <- unsafeRead iejeitjtisjsls 5
                              case (l /= is) && (l /= js) of
                                True -> do
                                  let it = minimum [l,is,js]
                                      jt = maximum [l,is,js]
                                  unsafeWrite iejeitjtisjsls 2 it
                                  unsafeWrite iejeitjtisjsls 3 jt
                                  let lt = is+js+l-it-jt
                                  dfr1 <- readArray frthdf (it,lt)
                                  dfr2 <- readArray frthdf (lt,jt)
                                  let dfr = dfr1 + dfr2
                                  case dfr >= x of
                                    True  -> loop (l+1) dfr l
                                    False -> loop (l+1) x ls
                                False -> loop (l+1) x ls
          !ls <- loop 1 0 0
          unsafeWrite iejeitjtisjsls 6 ls
          is <- unsafeRead iejeitjtisjsls 4
          js <- unsafeRead iejeitjtisjsls 5
          difil <- readArray frthdf (min is ls, max is ls)
          diflj <- readArray frthdf (min js ls, max js ls)
          let dfnx' = max difil diflj
          unsafeWrite dfmxdfnx 1 dfnx'
          case dfmx/cuttb < dfnx' && difil > diflj of
            True -> do
              is' <- unsafeRead iejeitjtisjsls 4
              js' <- unsafeRead iejeitjtisjsls 5
              it <- unsafeRead iejeitjtisjsls 2
              unsafeWrite iejeitjtisjsls 2 is'
              unsafeWrite iejeitjtisjsls 4 js'
              unsafeWrite iejeitjtisjsls 5 it
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
  is <- unsafeRead iejeitjtisjsls 4
  js <- unsafeRead iejeitjtisjsls 5
  vtiIO <- extractColumn v is
  vtjIO <- extractColumn v js
  vti <- array1dToUVectorD vtiIO
  vtj <- array1dToUVectorD vtjIO
  case nregions == 4 of
    True -> do
      let vt = UV.zipWith (((*0.5).).(+)) vti vtj
      replaceDimensions vrts2 (js,top) (js,sbs+1) vt
      replaceDimension vrts2 (is,sbs+1) vti
      replaceDimensions vrts2 (is,sbs+2) (is,sbs+3) vt
      replaceDimensions vrts2 (js,sbs+2) (js,sbs+3) vtj
      it <- unsafeRead iejeitjtisjsls 2
      jt <- unsafeRead iejeitjtisjsls 3
      (vtit, vtjt) <- getVectors nd vrts2 top it jt
      -- vtit <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,it,top)) vrts2)
      -- vtjt <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,jt,top)) vrts2)
      -- vtit <- array1dToUVectorD vtitIO -- (=<<) (return . fromList) (getElems vtitIO)
      -- vtjt <- array1dToUVectorD vtjtIO -- (=<<) (return . fromList) (getElems vtjtIO)
      let vtt = UV.zipWith (((*0.5).).(+)) vtit vtjt
      replaceDimensions vrts2 (jt,top) (it,sbs+1) vtt
      replaceDimension vrts2 (jt,sbs+1) vtjt
      (vti2, vtj2) <- getVectors nd vrts2 (sbs+2) it jt
      -- vti2 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,it,sbs+2)) vrts2)
      -- vtj2 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,jt,sbs+2)) vrts2) -- :: IO IO1dArray
      -- vti2 <- (=<<) (return . fromList) (getElems vti2IO)
      -- vtj2 <- (=<<) (return . fromList) (getElems vtj2IO)
      let vt2 = UV.zipWith (((*0.5).).(+)) vti2 vtj2
      replaceDimensions vrts2 (jt,sbs+2) (it,sbs+3) vt2
      replaceDimension vrts2 (jt,sbs+3) vtj2
    False -> do
      let vt = UV.zipWith (\x -> \y -> (2*x+y)/3) vti vtj
      replaceDimensions vrts2 (js,top) (is,sbs+1) vt
      case dfmx/cuttf < dfnx of
        True -> do
          replaceDimensions vrts2 (js,sbs+1) (is,sbs+2) vtj
          replaceDimension vrts2 (js,sbs+2) vtj
          ls <- unsafeRead iejeitjtisjsls 6
          (vtj1, vtl1) <- getVectors nd vrts2 (sbs+1) js ls
          -- vtj1 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,js,sbs+1)) vrts2)
          -- vtl1 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,ls,sbs+1)) vrts2)
          -- vtj1 <- (=<<) (return . fromList) (getElems vtj1IO)
          -- vtl1 <- (=<<) (return . fromList) (getElems vtl1IO)
          let vt1 = UV.zipWith (((*0.5).).(+)) vtj1 vtl1
          replaceDimensions vrts2 (ls,sbs+1) (js,sbs+2) vt1
          replaceDimension vrts2 (ls,sbs+2) vtl1
        False -> do
          let vv = UV.zipWith (\x -> \y -> (x+2*y)/3) vti vtj
          replaceDimensions vrts2 (js,sbs+1) (is,sbs+2) vv
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

replaceDimensions :: IO3dArray -> (Int,Int) -> (Int,Int) -> UVectorD -> IO ()
replaceDimensions m (j1,k1) (j2,k2) v = do
  (_, (n,_,_)) <- getBounds m
  let loop :: Int -> IO ()
      loop i | i == n+1 = return ()
             | otherwise = do
               writeArray m (i,j1,k1) ((UV.!) v (i-1))
               writeArray m (i,j2,k2) ((UV.!) v (i-1))
               loop (i+1)
  loop 1

-- | Number of evaluations for each subregion
smpchc :: Int -> Int -> Int
smpchc nd key =
  if key == 3
    then div ((nd+4)*(nd+3)*(nd+2)) 6 + (nd+2)*(nd+1)
    else if key == 1
      then 2*nd + 3
      else if key == 2
        then div ((nd+3)*(nd+2)) 2 + 2*(nd+1)
        else div ((nd+5)*(nd+4)*(nd+3)*(nd+2)) 24 + 5*(div ((nd+2)*(nd+1)) 2)

-- | Checks validity of parameters
check :: Int -> Int -> Int -> Double -> Double -> Int -> Int -> Int
check nd nf mxfs ea er sbs key =
  if ea < 0 || er < 0
    then 5
    else if nf < 1
      then 4
      else if nd < 2
        then 3
        else if key < 1 || key > 4
          then 2
          else if mxfs < sbs * smpchc nd key
            then 1
            else 0

adsimp :: Int -> Int -> Int -> (UVectorD -> UVectorD) -> Double -> Double
       -> Int -> IO3dArray -> Bool -> IO (UVectorD, UVectorD, Int, Bool)
adsimp nd nf mxfs f ea er key vrts info = do
  (_, (_,_,sbs)) <- getBounds vrts
  case check nd nf mxfs ea er sbs key of
    0 -> smpsad nd nf f mxfs ea er key (smpchc nd key) sbs vrts info
    1 -> smpsad nd nf f (sbs*(smpchc nd key)) ea er key (smpchc nd key) sbs vrts info
    2 -> error "integration rule must be between 1 and 4"
    3 -> error "dimension must be at least 2"
    4 -> error "number of components must be at least 1"
    5 -> error "requested errors must be positive"

type Params = (Bool, Int, Int, Seq UVectorD, Seq UVectorD,
               UVectorD, UVectorD, IO3dArray, Seq Double)

smpsad :: Int -> Int -> (UVectorD -> UVectorD) -> Int -> Double -> Double -> Int
       -> Int -> Int -> IO3dArray -> Bool -> IO (UVectorD, UVectorD, Int, Bool)
smpsad nd nf f mxfs ea er key rcls sbs vrts info = do
  let dfcost = 1 + 2*nd*(nd+1)
  (g, w, pospts) <- smprms nd key
  simplices <- mapM (toSimplex vrts (nd+1)) [1..sbs]
  let vol = S.fromList $ map simplexVolume simplices
      nv = sbs*rcls
  matrices <- mapM
              (\k -> mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,k)) vrts)
              (S.fromList [1..sbs])
  br <- mapM (\(m,v) -> smprul m nf f v g w pospts) (S.zip matrices vol)
  aes <- mapM (array1dToUVectorD.snd) br
  vls <- mapM (array1dToUVectorD.fst) br
  let vl = foldl1 (UV.zipWith (+)) vls
      ae = foldl1 (UV.zipWith (+)) aes
      fl = getFL ae vl
  let loop :: Params -> IO (UVectorD, UVectorD, Int, Bool)
      loop !params | not fl' || nv'+dfcost+4*rcls > mxfs =
                     return (vl', ae', nv', fl')
                   | otherwise = do
                      let maxs = fmap UV.maximum aes'
                          imax = fromJust $ S.findIndexL (== maximum maxs) maxs
                          vl0 = UV.zipWith (-) vl' (index vls' imax)
                          ae0 = UV.zipWith (-) ae' (index aes' imax)
                      (nregions, vrts2) <- smpdfs nd f (imax+1) sbs' vrts'
                      let vi = (index vol' imax)/(toDbl nregions)
                          nv2 = nv' + (nregions-1)*rcls -- nregions*rcls ?
                          sbs2 = sbs' + nregions-1
                      matrices2 <- mapM
                                   (\k -> mapIndices ((1,1),(nd,nd+1))
                                          (\(i,j) -> (i,j,k)) vrts2)
                                   (S.fromList ((imax+1):[(sbs'+1)..sbs2]))
                      br2 <- mapM (\m -> smprul m nf f vi g w pospts) matrices2
                      rgnerrs <- mapM (array1dToUVectorD.snd) br2
                      basvals <- mapM (array1dToUVectorD.fst) br2
                      let vl2 = UV.zipWith (+) vl0
                                (foldl1 (UV.zipWith (+)) basvals)
                          ae2 = UV.zipWith (+) ae0
                                (foldl1 (UV.zipWith (+)) rgnerrs)
                          aes2 = (update imax (index rgnerrs 0) aes')
                                 >< (S.drop 1 rgnerrs)
                          vls2 = (update imax (index basvals 0) vls')
                                 >< (S.drop 1 basvals)
                          fl2 = getFL ae2 vl2
                          vol2 = (update imax vi vol')
                                 >< (S.replicate (nregions-1) vi)
                      loop (fl2, nv2, sbs2, aes2, vls2, ae2, vl2, vrts2, vol2)
                    where
                      (fl',nv',sbs',aes',vls',ae',vl',vrts',vol') = params
  -- dans le code R il fait rowSums mais ça me semble inutile
  loop (fl, nv, sbs, aes, vls, ae, vl, vrts, vol)
  where
    getFL a v = UV.any (> (max ea (UV.maximum (UV.map ((*er).abs) v)))) a


toSimplex :: IO3dArray -> Int -> Int -> IO Simplex
toSimplex m n k = do
  (_, (nrow,_,_)) <- getBounds m
  let getColumn :: Int -> IO IO1dArray
      getColumn col = mapIndices (1,nrow) (\i -> (i,col,k)) m
  columns <- mapM getColumn [1..n]
  mapM getElems columns
