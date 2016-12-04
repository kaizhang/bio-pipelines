module Bio.Pipeline.ScanMotifs where

import           Bio.Data.Bed
import           Bio.Motif
import           Bio.Seq.IO
import           Conduit
import           Data.Default

scanMotifs :: FilePath   -- ^ Genome file
           -> [Motif]    -- ^ Motifs
           -> Double     -- ^ P-value
           -> FilePath   -- ^ Output
           -> [FilePath]
           -> IO FilePath
scanMotifs genome motifs p output fls = withGenome genome $ \g -> do
    beds <- mapM readBed' fls :: IO [[BED3]]
    mergeBed (concat beds) =$= motifScan g motifs def p =$=
        getMotifScore g motifs def =$= getMotifPValue (Just (1 - p * 10)) motifs def $$
        writeBed output
    return output
