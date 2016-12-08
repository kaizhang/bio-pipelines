{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}

module Bio.Pipeline.CallPeaks
    ( CallPeakOpts(..)
    , CallPeakOptSetter
    , tmpDir
    , qValue
    , gSize
    , pair
    , defaultCallPeakOpts
    , callPeaks
    ) where

import           Bio.Data.Experiment.Types
import           Control.Lens
import           Control.Monad.State.Lazy
import qualified Data.Text                 as T
import           Shelly                    hiding (FilePath)
import           System.IO.Temp            (withTempDirectory)

import           Bio.Pipeline.Utils

type CallPeakOptSetter = State CallPeakOpts ()

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir :: FilePath
    , callPeakOptsQValue :: Double
    , callPeakOptsGSize  :: String
    , callPeakOptsPair   :: Bool
    --, callPeakOptsBroad :: Bool
    --, callPeakOptsBroadCutoff :: Double
    }

makeFields ''CallPeakOpts

defaultCallPeakOpts :: CallPeakOpts
defaultCallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir = "./"
    , callPeakOptsQValue = 0.01
    , callPeakOptsGSize = "mm"
    , callPeakOptsPair  = False
    --, callPeakOptsBroad = False
    --, callPeakOptsBroadCutoff = 0.05
    }

-- | Call peaks using MACS2.
callPeaks :: FilePath             -- ^ Ouptut file
          -> [FileSet]            -- ^ one or more samples. Samples will be concatenated.
          -> [FileSet]            -- ^ zero or more input samples
          -> CallPeakOptSetter    -- ^ Options
          -> IO FileSet
callPeaks output targets' inputs' setter = do
    macs2 output (map (^.location) targets) (map (^.location) inputs)
        fileFormat opt
    return $ Single $ format .~ NarrowPeakFile $ location .~ output $
        tags .~ ["macs2"] $ emptyFile
  where
    targets = targets'^..folded._Single
    inputs = inputs'^..folded._Single
    opt = execState setter defaultCallPeakOpts
    fileFormat
        | opt^.pair = if allFilesEqual
            then case formt of
                BamFile -> "BAMPE"
                _ -> error "Only BAM input is supported in pairedend mode."
            else error "Only BAM input is supported in pairedend mode."
        | otherwise = "AUTO"
      where
        formt = head targets ^. format
        allFilesEqual = allEqual targets && allEqual inputs
          where
            allEqual = all (\x -> x^.format == formt)
{-# INLINE callPeaks #-}

macs2 :: FilePath      -- ^ Output
      -> [FilePath]    -- ^ Target
      -> [FilePath]    -- ^ Input
      -> String        -- ^ File format
      -> CallPeakOpts
      -> IO ()
macs2 output targets inputs fileformat opt = withTempDirectory (opt^.tmpDir)
    "tmp_macs2_dir." $ \tmp -> shelly $ silently $ do
        run_ "macs2" $
            [ "callpeak", "-f", T.pack fileformat, "-g", T.pack $ opt^.gSize
            , "--outdir", T.pack tmp, "--tempdir", T.pack tmp, "--keep-dup"
            , "all", "-q", T.pack $ show $ opt^.qValue
            ] ++ samples ++ controls
        mv (fromText $ T.pack $ tmp ++ "/NA_peaks.narrowPeak") $ fromText $
            T.pack output
  where
    controls | null inputs = []
             | otherwise = "-c" : map T.pack inputs
    samples | null targets = error "Empty sample."
            | otherwise = "-t" : map T.pack targets
{-# INLINE macs2 #-}
