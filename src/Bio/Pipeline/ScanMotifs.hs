module Bio.Pipeline.ScanMotifs where

import           Bio.Data.Bed
import           Bio.Motif
import           Bio.Seq.IO
import           Conduit
import           Data.Default

scanMotifs :: FilePath   -- ^ genome file
           -> FilePath   -- ^ motif file
           -> Double     -- ^ p-value
           -> FilePath   -- ^ output
           -> [FilePath]
           -> IO FilePath
scanMotifs genome motifs p output fls = withGenome genome $ \g -> do
    ms <- readMEME motifs
    beds <- mapM readBed' fls :: IO [[BED3]]
    sites <- mergeBed (concat beds) =$= motifScan g ms def p =$=
        getMotifScore g ms def $$ sinkList
    writeBed' output $ getMotifPValue ms def sites
    return output
