{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}
module Internal4
  where
import           Control.Monad               ((<$!>), (=<<))
import           Control.Monad.ST            (ST)
import           Data.Array.IArray           (Array)
import           Data.Array.ST               (STUArray, getBounds, getElems,
                                              mapIndices, newArray, newArray_,
                                              newListArray, readArray,
                                              writeArray)
import           Data.Array.Unboxed          (UArray)
import qualified Data.Array.Unboxed          as UA
import qualified Data.Foldable               as DF
import           Data.List                   (sort)
import           Data.Maybe                  (fromJust)
import           Data.Sequence               (Seq, index, update, (><))
import qualified Data.Sequence               as S
import qualified Data.Vector                 as V
import           Data.Vector.Unboxed         (Vector, empty, fromList, toList,
                                              unsafeFreeze)
import qualified Data.Vector.Unboxed         as UV
import           Data.Vector.Unboxed.Mutable (STVector, new, unsafeRead,
                                              unsafeWrite)
import qualified Data.Vector.Unboxed.Mutable as UMV
--import           Math.Combinat.Permutations  (permuteMultiset)
import           Simplex                     (Simplex, simplexVolume)

type STMatrix s = STUArray s (Int,Int) Double
type ST1dArray s = STUArray s Int Double
type ST3dArray s = STUArray s (Int,Int,Int) Double
type IMatrix = Array (Int,Int) Double
type UMatrix = UArray (Int,Int) Double
type U3dArray = UArray (Int,Int,Int) Double
type U1dArray = UArray Int Double
type STVectorD s = STVector s Double
type STVectorI s = STVector s Int
type UVectorD = Vector Double
type UVectorI = Vector Int

fromInt :: Int -> Double
fromInt = fromIntegral

smprms :: Int -> Int -> ST s (STMatrix s, Seq UVectorD, UVectorI)
smprms n key = do
  let (rls, gms, wts) | key == 1 = (3, 2, 3) :: (Int, Int, Int)
                      | key == 2 = (5, 4, 6) :: (Int, Int, Int)
                      | key == 3 = (7, 7, 11) :: (Int, Int, Int)
                      | key == 4 = (if n == 2 then (7, 11, 20) else (7, 12, 21))
                                    :: (Int, Int, Int)
  w <- newArray ((1,1), (wts, rls)) 0 :: ST s (STMatrix s)
  pts <- UMV.replicate wts 0 :: ST s (STVector s Int)
  g <- newArray ((1,1), (n+1, wts)) 0 :: ST s (STMatrix s)
  let np = n+1
      n2 = np * (n+2)
      n4 = n2 * (n+3) * (n+4)
      n6 = n4 * (n+5) * (n+6)
      n8 = n6 * (n+7) * (n+8)
      ndbl = fromInt n
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
  let n2double = fromInt n2 :: Double
  writeArray w (2,iw) ((ndbl+3)^2 / (4*n2double))
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
          unsafeWrite pts gms 3
          writeArray w (gms+1,iw-1) (1/6) -- rq: déjà écrit
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
      let n4double = fromInt n4
      writeArray w (2,iw-2) (-(ndbl+3)^5 / (16*n4double))
      writeArray w (3,iw-2) ((ndbl+5)^5 / (16*n4double*(ndbl+5)))
      writeArray w (4,iw-2) ((ndbl+5)^5 / (16*n4double*(ndbl+5)))
    False -> return ()
  case key > 2 of
    True -> do
      let u1 = (ndbl+7+2*sqrt15) / (ndbl^2+14*ndbl-11)
          v1 = (1-(ndbl-1)*u1)/2
          d1 = v1 - u1
      writeArray g (1,gms+3) v1
      writeArray g (2,gms+3) v1
      _ <- mapM (\j -> writeArray g (j,gms+3) u1) [3 .. np]
      unsafeWrite pts (gms+2) (div (n*np) 2)
      let u2 = (ndbl+7-2*sqrt15) / (ndbl^2+14*ndbl-11)
          v2 = (1-(ndbl-1)*u2)/2
          d2 = v2 - u2 -- utilisé? oui plus bas
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
              writeArray w (gms+3,iw-3) (2*15/567)
              unsafeWrite pts (gms+3) 0
            False -> do
              let r2 = (ndbl+4+sqrt15) / (ndbl*ndbl+8*ndbl+1)
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
          r = 2*sqrt(-p) -- 2*sqrt(-p*p*p)**(1/3)
          tp = 2*pi/3
          a1 = -a + r*(cos(th))
          a2 = -a + r*(cos(th+2*tp))
          a3 = -a + r*(cos(th+tp))
          npdbl = fromInt np
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
                 ((u7-(a2+a3)*u6+a2*a3*u5)/(a1^2-(a2+a3)*a1+a2*a3)/a1^5)
      writeArray w (gms+6,iw-5)
                 ((u7-(a1+a3)*u6+a1*a3*u5)/(a2^2-(a1+a3)*a2+a1*a3)/a2^5)
      writeArray w (gms+7,iw-5)
                 ((u7-(a2+a1)*u6+a2*a1*u5)/(a3^2-(a2+a1)*a3+a2*a1)/a3^5)
      writeArray g (1,gms+8) (4/(ndbl+7))
      writeArray g (2,gms+8) (4/(ndbl+7))
      _ <- mapM (\i -> writeArray g (i,gms+8) (1/(ndbl+7))) [3..np]
      unsafeWrite pts (gms+7) (div (np*n) 2)
      writeArray w (gms+8,iw-5) (10*(ndbl+7)^6/(fromInt $ 729*n6))
      writeArray g (1,gms+9) (11/(ndbl+7)/2)
      writeArray g (2,gms+9) (5/(ndbl+7)/2)
      _ <- mapM (\i -> writeArray g (i,gms+9) (1/(ndbl+7))) [3..np]
      unsafeWrite pts (gms+8) (np*n)
      writeArray w (gms+9,iw-5) (64*(ndbl+7)^6/(fromInt $ 6561*n6))
      writeArray w (4,iw-5) ((-(ndbl+5)^7)/(fromInt $ 64*n6))
      writeArray w (7,iw-5) (((ndbl+7)^7)/(fromInt $ 64*n6*(n+7)))
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
      writeArray w (2,iw-6) (- (ndbl+3)^9/(fromInt $ 6*256*n6))
      _ <- mapM (\i -> writeArray w (i,iw-6) ((ndbl+5)^9/(fromInt $ 512*n6*(n+7)))) [3,4]
      _ <- mapM (\i -> writeArray w (i,iw-6) (-(ndbl+7)^9/(fromInt $ 256*n8))) [5,6,7]
      _ <- mapM (\i -> writeArray w (i,iw-6) ((ndbl+9)^9/(fromInt $ 256*n8*(n+9)))) [8..11]
      case n > 2 of
        True ->  do
          _ <- mapM (\i -> writeArray g (i,12) (3/(ndbl+9))) [1..4]
          _ <- mapM (\i -> writeArray g (i,12) (1/(ndbl+9))) [5..np]
          unsafeWrite pts 11 (div (np*n*(n-1)*(n-2)) 24)
          writeArray w (12,iw-6) ((ndbl+9)^9/(fromInt $ 256*n8*(n+9)))
        False -> return ()
    False -> return ()
  -- uw <- unsafeFreeze w
  rowsIO <- mapM (extractRow w) (S.fromList [2..wts])
  rows <- mapM array1dToUVectorD rowsIO
  let cols = transpose rows
  pts_out <- unsafeFreeze pts
  let ptsU = UV.map fromInt pts_out
      row1 = fmap (\col -> 1 - UV.foldr (+) 0
                           (UV.zipWith (*) (UV.tail ptsU) col))
                  (S.take rls cols)
      wcols = fmap (\j -> UV.cons (index row1 j) (index cols j))
                   (S.fromList [0..(rls-1)])
      col1 = index wcols 0
      wcols2 = (S.<|) col1
                      (fmap (\col -> UV.zipWith (-) col col1) (S.drop 1 wcols))
      col2 = index wcols2 1
      nb = UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map (^2) col1))
      ratio = nb / (UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map (^2) col2)))
      wcol2 = UV.map (*(sqrt ratio)) col2
      wcols3 = (S.<|) col1 ((S.<|) wcol2 (S.drop 2 wcols2))
  let updateW :: Seq UVectorD -> Int -> Seq UVectorD
      updateW cols k = update (k-1) wknew cols
         where
          ptsW = UV.map (/nb) (UV.zipWith (*) ptsU (index cols (k-1)))
          slice = S.drop 1 (S.take (k-1) cols)
          prod1 = (fromList.(DF.toList)) $ -- ou alors pas de vecteurs
                  fmap ((UV.foldr (+) 0).(UV.zipWith (*) ptsW)) slice
          rows = transpose slice
          prod2 = (fromList.(DF.toList)) $
                  fmap ((UV.foldr (+) 0).(UV.zipWith (*) prod1)) rows
          wk = UV.zipWith (-) (index cols (k-1)) prod2
          ratio = nb / (UV.foldr (+) 0 (UV.zipWith (*) ptsU (UV.map (^2) wk)))
          wknew = UV.map (*(sqrt ratio)) wk
      wcolsnew = foldl updateW wcols3 [3..rls]
  return (g, transpose wcolsnew, pts_out)

transpose :: Seq UVectorD -> Seq UVectorD
transpose cols =
  fmap (\i -> (fromList.(DF.toList)) $ fmap (\col -> col UV.!i) cols)
       (S.fromList ([0..(UV.length (index cols 0) - 1)]))

matprod :: STMatrix s -> UVectorD -> ST s (UVectorD)
matprod mat x = do
  (_,(m,n)) <- getBounds mat
  out <- new m :: ST s (STVectorD s)
  let step i | i == m+1 = unsafeFreeze out
             | otherwise = do
              !coef <- innerstep i 1 0
              unsafeWrite out (i-1) coef
              step (i+1)
--      innerstep :: Int -> Int -> Double -> ST s (Double)
      innerstep i j !s | j == n+1 = return s
                       | otherwise = do
                        mat_ij <- readArray mat (i,j)
                        innerstep i (j+1) (s + mat_ij * (x UV.! (j-1)))
  step 1

-- matprod2 :: STMatrix s -> V.Vector UVectorD -> ST s ([UVectorD])
-- matprod2 mat x = do
--   (_,(m,n)) <- getBounds mat
--   let l = V.length x
--   out <- replicateM l (UMV.new m :: ST s (UMV.STVector s Double))
--   let step i | i == m+1 = mapM unsafeFreeze out
--              | otherwise = do
--               !coefs <- innerstep i 1 (UV.replicate l 0)
--               mapM (\k -> unsafeWrite (out!!k) (i-1) (coefs UV.! k)) [0..(l-1)]
--               step (i+1)
-- --      innerstep :: Int -> Int -> Double -> ST s (Double)
--       innerstep i j !s | j == n+1 = return s
--                        | otherwise = do
--                         mat_ij <- readArray mat (i,j)
--                         innerstep i (j+1)
--                                   (UV.zipWith (+) s (UV.map (*mat_ij) (UV.map (\k -> (x V.! k) UV.! (j-1)) (UV.fromList [0..(l-1)]))))
-- --                        innerstep i (j+1) (s + mat_ij * (x UV.! (j-1)))
--   step 1

unfold1 :: (UVectorD -> Maybe UVectorD) -> UVectorD -> V.Vector UVectorD
unfold1 f x = case f x of
  Nothing -> V.singleton x
  Just y  -> V.cons x (unfold1 f y)

-- permuteV :: [Double] -> [UVectorD]
-- permuteV xs = unfold1 next (fromList (sort xs)) where
--   next xs = case findj (UV.reverse xs, UV.empty) of
--     Nothing -> Nothing
--     Just (lls, rs) -> Just $ inc (UV.head lls) (UV.tail lls) (UV.reverse rs, UV.empty)
--   findj (xxs , yys) | xxs == empty = Nothing
--                     | yys == empty = findj (UV.tail xxs , UV.take 1 xxs)
--                     | otherwise = let x = UV.head xxs in if x >= (UV.head yys)
--                        then findj (UV.tail xxs , UV.cons x  yys )
--                        else Just (xxs, yys)
--   inc !u us (xxs , yys) = let x = UV.head xxs in if u >= x
--     then inc u us (UV.tail xxs , UV.cons x yys)
--     else UV.reverse (UV.cons x us) UV.++ UV.reverse (UV.cons u yys) UV.++ (UV.tail xxs)
--
permuteV2 :: [Double] -> V.Vector UVectorD
permuteV2 xs = unfold1 next (fromList (sort xs)) where
  next xs = case findj (UV.reverse xs, UV.empty) of
    Nothing -> Nothing
    Just (lls, rs) -> Just $ inc (UV.head lls) (UV.tail lls) (UV.reverse rs, UV.empty)
  findj (xxs , yys) | xxs == empty = Nothing
                    | yys == empty = findj (UV.tail xxs , UV.take 1 xxs)
                    | otherwise = let x = UV.head xxs in if x >= (UV.head yys)
                       then findj (UV.tail xxs , UV.cons x  yys )
                       else Just (xxs, yys)
  inc !u us (xxs , yys) = let x = UV.head xxs in if u >= x
    then inc u us (UV.tail xxs , UV.cons x yys)
    else UV.reverse (UV.cons x us) UV.++ UV.reverse (UV.cons u yys) UV.++ (UV.tail xxs)

smpsms :: STMatrix s -> Int -> (UVectorD -> UVectorD) -> ST1dArray s
       -> Double -> ST s (ST1dArray s)
smpsms vertex nf f g scalar = do
  gAsList <- getElems g
  -- (_,(nrow,ncol)) <- getBounds vertex -- on dirait que ce truc a ralenti
  f_gPermuts <- V.mapM ((fmap f).(matprod vertex))
                       (permuteV2 gAsList) -- c'est ça qui est time-consuming
  -- prods <- matprod2 vertex (permuteV2 gAsList)
  -- let f_gPermuts = map f prods
  newListArray (1,nf)
               (toList (UV.map (*scalar) (foldl1 (UV.zipWith (+)) f_gPermuts)))

extractColumn :: STMatrix s -> Int -> ST s (ST1dArray s)
extractColumn m j = do
 (_, (nrow,_)) <- getBounds m
 mapIndices (1,nrow) (\i -> (i,j)) m

extractRow :: STMatrix s -> Int -> ST s (ST1dArray s)
extractRow m i = do
 (_, (_,ncol)) <- getBounds m
 mapIndices (1,ncol) (\j -> (i,j)) m

outerProduct2 :: ST1dArray s -> UVectorD -> ST s (STMatrix s)
outerProduct2 x1 x2 = do
 (_, n1) <- getBounds x1
 let n2 = UV.length x2
 out <- newArray_ ((1,1),(n1,n2)) :: ST s (STMatrix s)
 --let step :: Int -> ST s (STMatrix s)
 let step i | i == n1+1 = return out
            | otherwise = do
               x1_i <- readArray x1 i
               inner x1_i 0
             where
               inner !x j | j == n2 = step (i+1)
                          | otherwise = do
                             writeArray out (i,j+1) (x * (x2 UV.! j))
                             inner x (j+1)
 step 1

sumMatrices :: [STMatrix s] -> ST s (STMatrix s)
sumMatrices matrices = do
 (_, (n1,n2)) <- getBounds (head matrices)
 out <- newArray_ ((1,1),(n1,n2)) :: ST s (STMatrix s)
 let step i | i == n1+1 = return out
            | otherwise = inner 1
              where
                inner j | j == n2+1 = step (i+1)
                        | otherwise = do
                          coefs <- mapM (\m -> readArray m (i,j)) matrices
                          writeArray out (i,j) (sum coefs)
                          inner (j+1)
 -- let step i j | i == n1 && j == n2+1 = return out
 --              | j == n2+1 = step (i+1) 1 -- avec un inner tu gagnerais un test
 --              | otherwise = do
 --                coefs <- mapM (\m -> readArray m (i,j)) matrices
 --                writeArray out (i,j) (sum coefs)
 --                step i (j+1)
 step 1

smprul :: STMatrix s -> Int -> (UVectorD -> UVectorD) -> Double -> STMatrix s
      -> Seq UVectorD -> UVectorI -> ST s (ST1dArray s, ST1dArray s)
smprul vrts nf f vol g w pts = do
 let rtmn = 0.1
     small = 1e-12
     errcof = 8
     rls = UV.length (index w 0)
     ptsPositive = toList $ UV.findIndices (> 0) pts
 toSum <- mapM (\k -> do
                        g_colk <- extractColumn g (k+1)
                        sms <- smpsms vrts nf f g_colk vol
                        outerProduct2 sms (index w k))
               ptsPositive
 rule <- sumMatrices toSum -- quid si la liste est vide ?
 basval <- extractColumn rule 1
 rgnerr <- newArray (1,nf) 0 :: ST s (ST1dArray s)
 --let step :: Int -> ST s ()
 let step i | i == nf+1 = return ()
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
              --innerstep :: Int -> Double -> Double -> ST s (Double, Double)
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

rowMeans :: STMatrix s -> ST s UVectorD
rowMeans m = do
 (_, (nrow,ncol)) <- getBounds m
 outIO <- new nrow :: ST s (STVectorD s)
 --let step :: Int -> ST s ()
 let step i | i == nrow+1 = return ()
            | otherwise = do
               !sum_i <- inner 1 0
               unsafeWrite outIO (i-1) sum_i
               step (i+1)
             where
               --inner :: Int -> Double -> ST s Double
               inner j !x | j == ncol+1 = return (x / (fromInt ncol))
                          | otherwise = do
                            coef <- readArray m (i,j)
                            inner (j+1) (x + coef)
 step 1
 unsafeFreeze outIO

array1dToUVectorD :: ST1dArray s -> ST s UVectorD
array1dToUVectorD array = do
 (<$!>) fromList (getElems array)

smpdfs0 :: Int -> (UVectorD -> UVectorD) -> Int
       -> ST3dArray s -> ST s (STVectorI s, STVectorD s, Int)
smpdfs0 nd f top vrts = do -- nf et sbs pas utilisés
  let cuttf = 2.0
      cuttb = 8.0
  v <- mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,top)) vrts
  cn <- rowMeans v
  let fc = f cn
      dfmd = UV.foldr (+) 0 (UV.map abs fc)
  frthdf <- newArray ((1,1),(nd,nd+1)) 0 :: ST s (STMatrix s)
  iejeitjtisjsls <- new 7 :: ST s (STVectorI s)
  unsafeWrite iejeitjtisjsls 4 1
  unsafeWrite iejeitjtisjsls 5 2
  dfmxdfnx <- UMV.replicate 2 0 :: ST s (STVectorD s)
  --let step :: Int -> Double -> ST s ()
  let step i x | i == nd+1 = return ()
               | otherwise = do
                  !emx <- inner (i+1) x
                  step (i+1) emx
                where
--                  inner :: Int -> Double -> ST s (Double)
                   inner j !y | j == nd+2 = return y
                              | otherwise = do
                                vi <- (=<<) array1dToUVectorD (extractColumn v i)
                                vj <- (=<<) array1dToUVectorD (extractColumn v j)
                                let h = UV.map (*(2/(5*(fromInt nd +1))))
                                               (UV.zipWith (-) vi vj)
                                    ewd = UV.foldr (+) 0 (UV.map abs h)
                                    twoh = UV.map (*2) h
                                    t1 = f (UV.zipWith (-) cn twoh)
                                    t2 = f (UV.zipWith (+) cn twoh)
                                    t3 = UV.map (*6) fc
                                    t4 = f (UV.zipWith (-) cn h)
                                    t5 = f (UV.zipWith (+) cn h)
                                    t6 = UV.zipWith (((*(-4)).).(+)) t4 t5
                                    tsum = foldl1 (UV.zipWith (+)) [t1,t2,t3,t6]
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
          --let -- loop :: Int -> Double -> Int -> ST s (Int)
          let loop l !x !ls | l == nd+2 = return ls
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
          let dfnx = max difil diflj
          unsafeWrite dfmxdfnx 1 dfnx
          case dfmx/cuttb < dfnx && difil > diflj of
            True -> do
              is <- unsafeRead iejeitjtisjsls 4
              js <- unsafeRead iejeitjtisjsls 5
              it <- unsafeRead iejeitjtisjsls 2
              unsafeWrite iejeitjtisjsls 2 is
              unsafeWrite iejeitjtisjsls 4 js
              unsafeWrite iejeitjtisjsls 5 it
            False -> return ()
  return (iejeitjtisjsls, dfmxdfnx, nregions)

smpdfs :: Int -> Int -> Int
       -> ST3dArray s -> (STVectorI s, STVectorD s, Int) -> ST s (ST3dArray s)
smpdfs nd top sbs vrts xxx = do
  let cuttf = 2.0
  v <- mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,top)) vrts
  let (iejeitjtisjsls, dfmxdfnx, nregions) = xxx
  vrts2 <- newArray_ ((1,1,1),(nd,nd+1,sbs+nregions-1)) :: ST s (ST3dArray s)
  let go1 i | i == nd+1 = return ()
            | otherwise = do
              --let go2 :: Int -> ST s ()
              let go2 j | j == nd+2 = go1 (i+1)
                        | otherwise = do
                          -- let go3 :: Int -> ST s ()
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
  vti <- (=<<) array1dToUVectorD (extractColumn v is)
  vtj <- (=<<) array1dToUVectorD (extractColumn v js)
  case nregions == 4 of
    True -> do
      let vt = UV.zipWith (((*0.5).).(+)) vti vtj
      replaceDimension vrts2 (js,top) vt
      replaceDimension vrts2 (is,sbs+1) vti
      replaceDimension vrts2 (js,sbs+1) vt
      replaceDimension vrts2 (is,sbs+2) vt
      replaceDimension vrts2 (is,sbs+3) vt
      replaceDimension vrts2 (js,sbs+2) vtj
      replaceDimension vrts2 (js,sbs+3) vtj
      it <- unsafeRead iejeitjtisjsls 2
      jt <- unsafeRead iejeitjtisjsls 3
      vtit <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,it,top)) vrts2)
      vtjt <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,jt,top)) vrts2)
      let vtt = UV.zipWith (((*0.5).).(+)) vtit vtjt
      replaceDimension vrts2 (jt,top) vtt
      replaceDimension vrts2 (it,sbs+1) vtt
      replaceDimension vrts2 (jt,sbs+1) vtjt
      vti2 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,it,sbs+2)) vrts2) -- :: IO IO1dArray
      vtj2 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,jt,sbs+2)) vrts2)  -- :: IO IO1dArray
      let vt2 = UV.zipWith (((*0.5).).(+)) vti2 vtj2
      replaceDimension vrts2 (jt,sbs+2) vt2
      replaceDimension vrts2 (it,sbs+3) vt2
      replaceDimension vrts2 (jt,sbs+3) vtj2
    False -> do
      let vt = UV.zipWith (\x -> \y -> (2*x+y)/3) vti vtj
      replaceDimension vrts2 (js,top) vt
      replaceDimension vrts2 (is,sbs+1) vt
      dfmx <- unsafeRead dfmxdfnx 0
      dfnx <- unsafeRead dfmxdfnx 1
      case dfmx/cuttf < dfnx of
        True -> do
          replaceDimension vrts2 (js,sbs+1) vtj
          replaceDimension vrts2 (is,sbs+2) vt
          replaceDimension vrts2 (js,sbs+2) vtj
          ls <- unsafeRead iejeitjtisjsls 6
          vtj1 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,js,sbs+1)) vrts2)
          vtl1 <- (=<<) array1dToUVectorD (mapIndices (1,nd) (\i -> (i,ls,sbs+1)) vrts2)
          let vt1 = UV.zipWith (((*0.5).).(+)) vtj1 vtl1
          replaceDimension vrts2 (ls,sbs+1) vt1
          replaceDimension vrts2 (js,sbs+2) vt1
          replaceDimension vrts2 (ls,sbs+2) vtl1
        False -> do
          let vv = UV.zipWith (\x -> \y -> (x+2*y)/3) vti vtj
          replaceDimension vrts2 (js,sbs+1) vv
          replaceDimension vrts2 (is,sbs+2) vv
          replaceDimension vrts2 (js,sbs+2) vtj
--   vti <- extractColumn v is
--   vtj <- extractColumn v js
--   case nregions == 4 of
--     True -> do
--       vt <- zip1dArraysWith vti vtj (\x -> \y -> (x+y)/2)
--       replaceDimensionIO vrts2 (js,top) vt
--       replaceDimensionU vrts2 (is,sbs+1) vti
--       replaceDimensionIO vrts2 (js,sbs+1) vt
--       replaceDimensionIO vrts2 (is,sbs+2) vt
--       replaceDimensionIO vrts2 (is,sbs+3) vt
--       replaceDimensionU vrts2 (js,sbs+2) vtj
--       replaceDimensionU vrts2 (js,sbs+3) vtj
--       it <- unsafeRead iejeitjtisjsls 2
--       jt <- unsafeRead iejeitjtisjsls 3
-- --      vtit <- mapIndices (1,nd) (\i -> (i,it,top)) vrts2 -- :: IO IO1dArray
--       vtjt <- mapIndices (1,nd) (\i -> (i,jt,top)) vrts2 -- :: IO IO1dArray
-- --      let vtt = UV.map (/2) (UV.zipWith (+) vtit vtjt)
--       -- vtit_f <- STA.freeze vtit
--       -- vtjt_f <- STA.freeze vtjt
--       -- vtt <- zip1dArraysWith vtit_f vtjt_f (\x -> \y -> (x+y)/2)
--       --replaceDimensionIO vrts2 (jt,top) vtt
--       replaceDimensionWith vrts2 (jt,top) (it,sbs+1) it jt top top (\x -> \y -> (x+y)/2)
-- --      replaceDimensionIO vrts2 (it,sbs+1) vtt
--       replaceDimensionIO vrts2 (jt,sbs+1) vtjt
--       -- vti2IO <- mapIndices (1,nd) (\i -> (i,it,sbs+2)) vrts2 -- :: IO IO1dArray
--       vtj2IO <- mapIndices (1,nd) (\i -> (i,jt,sbs+2)) vrts2  -- :: IO IO1dArray
--       -- vti2 <- (=<<) (return . fromList) (getElems vti2IO)
--       -- vtj2 <- (=<<) (return . fromList) (getElems vtj2IO)
-- --      let vt2 = UV.map (/2) (UV.zipWith (+) vti2 vtj2)
--       -- vti2IO_f <- STA.freeze vti2IO
--       -- vtj2IO_f <- STA.freeze vtj2IO
--       -- vt2 <- zip1dArraysWith vti2IO_f vtj2IO_f (\x -> \y -> (x+y)/2)
--       -- replaceDimensionIO vrts2 (jt,sbs+2) vt2
--       -- replaceDimensionIO vrts2 (it,sbs+3) vt2
--       replaceDimensionWith vrts2 (jt,sbs+2) (it,sbs+3) it jt (sbs+2) (sbs+2) (\x -> \y -> (x+y)/2)
--       replaceDimensionIO vrts2 (jt,sbs+3) vtj2IO
-- --      replaceDimension'' vrts2 (jt,sbs+3) (sbs+2)
--     False -> do
-- --      let vt = UV.map (/3) (UV.zipWith (+) (UV.map (*2) vti) vtj)
--       vt <- zip1dArraysWith vti vtj (\x -> \y -> (2*x+y)/3)
--       replaceDimensionIO vrts2 (js,top) vt
--       replaceDimensionIO vrts2 (is,sbs+1) vt
--       dfmx <- unsafeRead dfmxdfnx 0
--       dfnx <- unsafeRead dfmxdfnx 1 -- tu peux freezer dans smpdfs0
--       case dfmx/cuttf < dfnx of
--         True -> do
--           replaceDimensionU vrts2 (js,sbs+1) vtj
--           replaceDimensionIO vrts2 (is,sbs+2) vt
--           replaceDimensionU vrts2 (js,sbs+2) vtj
--           ls <- unsafeRead iejeitjtisjsls 6
--           -- vtj1IO <- mapIndices (1,nd) (\i -> (i,js,sbs+1)) vrts2 :: IO IO1dArray
--           -- vtl1IO <- mapIndices (1,nd) (\i -> (i,ls,sbs+1)) vrts2 :: IO IO1dArray
-- --          vtj1 <- STA.mapIndices (1,nd) (\i -> (i,js,sbs+1)) vrts2 -- :: ST s2 (STA.STUArray s1 Int Double)
--           vtl1 <- STA.mapIndices (1,nd) (\i -> (i,ls,sbs+1)) vrts2 -- :: STA.STUArray s Int Double
-- --          let vt1 = UV.map (/2) (UV.zipWith (+) vtj1 vtl1)
--           -- vtj1_f <- STA.freeze vtj1
--           -- vtl1_f <- STA.freeze vtl1
--           -- vt1 <- zip1dArraysWith vtj1_f vtl1_f (\x -> \y -> (x+y)/2)
--           -- replaceDimensionIO vrts2 (ls,sbs+1) vt1
--           -- replaceDimensionIO vrts2 (js,sbs+2) vt1
--           replaceDimensionWith vrts2 (ls,sbs+1) (js,sbs+2) js ls (sbs+1) (sbs+1) (\x -> \y -> (x+y)/2)
--           replaceDimensionIO vrts2 (ls,sbs+2) vtl1
--         False -> do
-- --          let vv = UV.map (/3) (UV.zipWith (+) vti (UV.map (*2) vtj))
--           vv <- zip1dArraysWith vti vtj (\x -> \y -> (x+2*y)/3)
--           replaceDimensionIO vrts2 (js,sbs+1) vv
--           replaceDimensionIO vrts2 (is,sbs+2) vv
--           replaceDimensionU vrts2 (js,sbs+2) vtj
  return vrts2
  -- where
  --       -- replaceDimensionIO :: (STA.STUArray s (Int,Int,Int) Double) -> (Int,Int) -> STA.STUArray s Int Double -> ST s ()
  --       replaceDimensionWith m (j,k) (jj,kk) j1 j2 k1 k2 f = loop 1
  --                                         where
  --                                           loop i | i == nd+1 = return ()
  --                                                  | otherwise = do
  --                                                    coef1 <- STA.readArray m (i,j1,k1)
  --                                                    coef2 <- STA.readArray m (i,j2,k2)
  --                                                    STA.writeArray m (i,j,k) (f coef1 coef2)
  --                                                    STA.writeArray m (i,jj,kk) (f coef1 coef2)
  --                                                    loop (i+1)
  --       replaceDimension'' m (j,k) kk = loop 1
  --                                       where
  --                                         --loop :: Int -> ST s ()
  --                                         loop i | i == nd+1 = return ()
  --                                                | otherwise = do
  --                                                  coef <- STA.readArray m (i,j,kk)
  --                                                  STA.writeArray m (i,j,k) coef
  --                                                  loop (i+1)
  --       replaceDimensionIO m (j,k) v = loop 1
  --                                       where
  --                                         --loop :: Int -> ST s ()
  --                                         loop i | i == nd+1 = return ()
  --                                                | otherwise = do
  --                                                  coef <- STA.readArray v i
  --                                                  STA.writeArray m (i,j,k) coef
  --                                                  loop (i+1)
  --       -- replaceDimensionU :: STA.STUArray s (Int,Int,Int) Double -> (Int,Int) -> U1dArray -> ST s ()
  --       replaceDimensionU m (j,k) v = loop 1
  --            where
  --                -- loop :: Int -> ST s ()
  --                loop i | i == nd+1 = return ()
  --                       | otherwise = do
  --                         coef <- readArray v i
  --                         STA.writeArray m (i,j,k) coef
  --                         loop (i+1)

replaceDimension :: ST3dArray s -> (Int,Int) -> UVectorD -> ST s ()
replaceDimension m (j,k) v = do
  (_, (n,_,_)) <- getBounds m
--  let loop :: Int -> IO ()
  let loop i | i == n+1 = return ()
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

adsimp :: Int -> Int -> Int -> (UVectorD -> UVectorD) -> Double -> Double
       -> Int -> ST3dArray s -> Bool -> ST s (UVectorD, UVectorD, Int, Bool)
adsimp nd nf mxfs f ea er key vrts info = do
  case key == 0 of
    True -> adsimp nd nf mxfs f ea er 3 vrts info
    False -> do
      (_, (_,_,sbs)) <- getBounds vrts
      let b = smpchc nd key
      smpsad nd nf f mxfs ea er key b sbs vrts info

type Params s = (Bool, Int, Int, Seq UVectorD, Seq UVectorD,
               UVectorD, UVectorD, ST3dArray s, Seq Double)

smpsad :: Int -> Int -> (UVectorD -> UVectorD) -> Int -> Double -> Double -> Int
       -> Int -> Int -> ST3dArray s -> Bool -> ST s (UVectorD, UVectorD, Int, Bool)
smpsad nd nf f mxfs ea er key rcls sbs vrts info = do
  let dfcost = 1 + 2*nd*(nd+1)
  (g, w, pts) <- smprms nd key
  simplices <- mapM (toSimplex vrts (nd+1)) [1..sbs]
  let vol = S.fromList $ map simplexVolume simplices
      nv = sbs*rcls
  matrices <- mapM
              (\k -> mapIndices ((1,1),(nd,nd+1)) (\(i,j) -> (i,j,k)) vrts)
              (S.fromList [1..sbs])
  br <- mapM (\(m,v) -> smprul m nf f v g w pts) (S.zip matrices vol)
  aes <- mapM (array1dToUVectorD.snd) br
  vls <- mapM (array1dToUVectorD.fst) br
  let vl = foldl1 (UV.zipWith (+)) vls
      ae = foldl1 (UV.zipWith (+)) aes
      fl = UV.any (> (max ea (UV.maximum (UV.map ((*er).abs) vl)))) ae
  let loop !params | not fl || nv+dfcost+4*rcls > mxfs = return (vl, ae, nv, fl)
                   | otherwise = do
                      let maxs = fmap UV.maximum aes
                          id = fromJust $ S.findIndexL (== maximum maxs) maxs
                          vl0 = UV.zipWith (-) vl (index vls id)
                          ae0 = UV.zipWith (-) ae (index aes id)
                      xxx <- smpdfs0 nd f (id+1) vrts
                      let (_,_,nregions) = xxx
                      vrts2 <- smpdfs nd (id+1) sbs vrts xxx
                      let vi = (index vol id)/(fromInt nregions)
                          nv2 = nv + (nregions-1)*rcls -- nregions*rcls ?
                          sbs2 = sbs + nregions-1
                      matrices <- mapM
                                  (\k -> mapIndices ((1,1),(nd,nd+1))
                                         (\(i,j) -> (i,j,k)) vrts2)
                                  (S.fromList ((id+1):[(sbs+1)..sbs2]))
                      br <- mapM (\m -> smprul m nf f vi g w pts) matrices
                      rgnerrs' <- mapM (array1dToUVectorD.snd) br
                      basvals' <- mapM (array1dToUVectorD.fst) br
                      let vl2 = UV.zipWith (+) vl0
                                               (foldl1 (UV.zipWith (+)) basvals')
                          ae2 = UV.zipWith (+) ae0
                                               (foldl1 (UV.zipWith (+)) rgnerrs')
                          aes2 = (update id (index rgnerrs' 0) aes)
                                 >< (S.drop 1 rgnerrs')
                          vls2 = (update id (index basvals' 0) vls)
                                 >< (S.drop 1 basvals')
                          fl2 = UV.any
                                (> (max ea (UV.maximum (UV.map ((*er).abs) vl2))))
                                ae2
                          vol2 = (update id vi vol)
                                 >< (S.replicate (nregions-1) vi)
                      loop (fl2, nv2, sbs2, aes2, vls2, ae2, vl2, vrts2, vol2)
                    where
                      (fl, nv, sbs, aes, vls, ae, vl, vrts, vol) = params
  -- dans le code R il fait rowSums mais ça me semble inutile
  loop (fl, nv, sbs, aes, vls, ae, vl, vrts, vol)

--
zip1dArraysWith :: ST1dArray s -> ST1dArray s -> (Double -> Double -> Double) -> ST s (ST1dArray s)
zip1dArraysWith a1 a2 f = do
  (_, n) <- getBounds a1
  out <- newArray_ (1,n) :: ST s (ST1dArray s)
  loop 1 out n
  where
      --loop :: Int -> STA.STUArray s Int Double -> Int -> ST s (STA.STUArray s Int Double)
      loop i out n | i == n+1 = return out
                   | otherwise = do
                     coef1 <- readArray a1 i
                     coef2 <- readArray a2 i
                     writeArray out i (f coef1 coef2)
                     loop (i+1) out n

zip1dArraysWith' :: U1dArray -> U1dArray -> (Double -> Double -> Double) -> ST s (ST1dArray s)
zip1dArraysWith' a1 a2 f = do
  let (_, n) = UA.bounds a1
  out <- newArray_ (1,n) :: ST s (ST1dArray s)
  loop 1 out n
  where
      loop :: Int -> ST1dArray s -> Int -> ST s (ST1dArray s)
      loop i out n | i == n+1 = return out
                   | otherwise = do
                     writeArray out i (f (a1 UA.! i) (a2 UA.! i))
                     loop (i+1) out n
--
-- zip1dArraysWith' :: (forall s.(STA.STUArray s (Int,Int,Int) Double)) -> Int -> Int -> Int -> (Double -> Double -> Double) -> UArray Int Double
-- zip1dArraysWith' m n j k f = STA.runSTUArray $ do
--   --(_, n) <- getBounds a1
--   out <- newArray_ (1,n) -- :: ST s (STA.STUArray s Int Double)
--   -- a1 <- STA.mapIndices (1,n) (\i -> (i,j,k)) m
--   -- a2 <- STA.mapIndices (1,n) (\i -> (i,j,k)) m
--   loop 1 out n
--   where
--       loop :: Int -> STA.STUArray s Int Double -> Int -> (forall s1 s2.(ST s2 (STA.STUArray s1 Int Double)))
--       loop i out n | i == n+1 = return out
--                    | otherwise = do
--                      coef1 <- STA.readArray m (i,j,k)
--                      coef2 <- STA.readArray m (i,j,k)
--                      STA.writeArray out i (f coef1 coef2)
--                      loop (i+1) out n

toSimplex :: ST3dArray s -> Int -> Int -> ST s Simplex
toSimplex m n k = do
  (_, (nrow,_,_)) <- getBounds m
  --let getColumn :: Int -> ST s (ST1dArray s)
  let getColumn col = mapIndices (1,nrow) (\i -> (i,col,k)) m
  columns <- mapM getColumn [1..n]
  mapM getElems columns
