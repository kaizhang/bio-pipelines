{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}
{-# LANGUAGE TypeSynonymInstances   #-}

module Bio.Pipeline.CallPeaks where

import           Bio.Data.Experiment.Types
import           Control.Lens
import           Control.Monad.State.Lazy
import qualified Data.Text                 as T
import           Turtle                    hiding (FilePath, format)
import qualified Turtle                    as T

import Bio.Pipeline.Utils

type CallPeakOptSetter = State CallPeakOpts ()

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir :: FilePath
    , callPeakOptsQValue :: Double
    , callPeakOptsGSize  :: String
    --, callPeakOptsBroad :: Bool
    --, callPeakOptsBroadCutoff :: Double
    }

makeFields ''CallPeakOpts

defaultCallPeakOpts :: CallPeakOpts
defaultCallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir = "./"
    , callPeakOptsQValue = 0.01
    , callPeakOptsGSize = "mm"
    --, callPeakOptsBroad = False
    --, callPeakOptsBroadCutoff = 0.05
    }

-- | Input: A list of target-input duos. The first BED file in Experiment will be
-- analyzed.
callPeaks :: FilePath
          -> CallPeakOptSetter
          -> [(Experiment, Maybe Experiment)]
          -> IO [Experiment]
callPeaks dir setter datasets = forM datasets $ \(target, input) -> do
    let targetFile = getFile target
        inputFile = fmap getFile input
        output = T.unpack $ T.format (s%"/"%s%"_rep"%d%".NarrowPeak") (T.pack dir)
            (target^.eid) (targetFile^.replication)
        peakFile = format .~ NarrowPeakFile $
            location .~ output $
            keywords .~ ["macs2"] $ targetFile
    callPeaksHelper targetFile inputFile output opt
    return $ files .~ [peakFile] $ target
  where
    opt = execState setter defaultCallPeakOpts
    getFile x = head $ filter
        (\fl -> fl^.format == BedFile || fl^.format == BedGZip) $ x^.files


callPeaksHelper :: File         -- ^ target
          -> Maybe File   -- ^ input
          -> FilePath
          -> CallPeakOpts
          -> IO ()
callPeaksHelper target input output opt = with ( mktempdir (fromText $ T.pack tmp)
    "macs2_tmp_dir_XXXXXX" ) $ \tmpDir -> do
        let cmd1 = T.format ("macs2 callpeak -t "%s%" -f BED -g "%s%" --outdir "%fp%
                " --tempdir "%fp%" --keep-dup all -q "%f)
                (T.pack $ getFile target) (T.pack $ opt^.gSize)
                tmpDir tmpDir (opt^.qValue)
            cmd2 = case input of
                Nothing -> ""
                Just input' -> T.format ("-c "%s) $ T.pack $ getFile input'
        shells (T.unwords [cmd1, cmd2]) empty
        mv (fromText $ T.format (fp%"/NA_peaks.narrowPeak") tmpDir) $
            fromText $ T.pack output
  where
    getFile fl | fl^.format == BedFile || fl^.format == BedGZip = fl^.location
               | otherwise = error "Incorrect file type"
    tmp = opt^.tmpDir
