{-# LANGUAGE ForeignFunctionInterface #-}
module Export
  where
import qualified Cubature as CUB
import Cubature (Result(..))
-- import qualified Cubature2 as CUB2
-- import qualified Cubature3 as CUB3
import qualified Cubature4 as CUB4
import           Foreign
import           Foreign.C
import Internal (smprms)
import qualified Data.Foldable               as DF
import qualified Data.Vector.Unboxed         as UV

-- import Internal (seqToVector, seqToVector2)
-- import Data.Vector.Unboxed as UV
-- import Data.Sequence as S
--
-- foreign export ccall ss1 :: Ptr CInt -> Ptr CInt -> IO ()
-- ss1 :: Ptr CInt -> Ptr CInt -> IO ()
-- ss1 n result = do
--   n <- peek n
--   v <- seqToVector (S.fromList [1..(fromIntegral n)])
--   poke result $ fromIntegral (UV.length v)
-- --
-- foreign export ccall ss2 :: Ptr CInt -> Ptr CInt -> IO ()
-- ss2 :: Ptr CInt -> Ptr CInt -> IO ()
-- ss2 n result = do
--   n <- peek n
--   let v = seqToVector2 (S.fromList [1..(fromIntegral n)])
--   poke result $ fromIntegral (UV.length v)

foreign export ccall testw :: Ptr Double -> IO ()
testw :: Ptr Double -> IO ()
testw result = do
  (_,w,_) <- smprms 4 3
  pokeArray result $ concat $ map UV.toList (DF.toList w)

foreign export ccall test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test rule maxevals result = do
  rule <- peek rule
  maxevals <- peek maxevals
  r <- CUB.example2 (fromIntegral maxevals) (fromIntegral rule)
  pokeArray result [realToFrac (value r), realToFrac (errorEstimate r),
                    fromIntegral (evaluations (r::Result)),
                    (fromIntegral.fromEnum) (success (r::Result))]

-- foreign export ccall test2 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
-- test2 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
-- test2 rule maxevals result = do
--   rule <- peek rule
--   maxevals <- peek maxevals
--   ([value], [errest], nevals, fl)
--     <- CUB2.example2 (fromIntegral maxevals) (fromIntegral rule)
--   pokeArray result [realToFrac value, realToFrac errest,
--                     fromIntegral nevals, (fromIntegral.fromEnum) fl]
--
-- foreign export ccall test2p :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
-- test2p :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
-- test2p rule maxevals result = do
--   rule <- peek rule
--   maxevals <- peek maxevals
--   ([value], [errest], nevals, fl)
--     <- CUB2.example2' (fromIntegral maxevals) (fromIntegral rule)
--   pokeArray result [realToFrac value, realToFrac errest,
--                     fromIntegral nevals, (fromIntegral.fromEnum) fl]
--
-- foreign export ccall test3 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
-- test3 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
-- test3 rule maxevals result = do
--   rule <- peek rule
--   maxevals <- peek maxevals
--   let ([value], [errest], nevals, fl) = runST $ CUB3.example2 (fromIntegral maxevals) (fromIntegral rule)
--   pokeArray result [realToFrac value, realToFrac errest,
--                     fromIntegral nevals, (fromIntegral.fromEnum) fl]
--
foreign export ccall test4 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test4 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test4 rule maxevals result = do
  rule <- peek rule
  maxevals <- peek maxevals
  let r = CUB4.example2 (fromIntegral maxevals) (fromIntegral rule)
  pokeArray result [realToFrac (CUB4.value r), realToFrac (CUB4.errorEstimate r),
                    fromIntegral (CUB4.evaluations (r::CUB4.Result)),
                    (fromIntegral.fromEnum) (CUB4.success (r::CUB4.Result))]
